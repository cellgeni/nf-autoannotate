#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MARKER_FILE="${REPO_ROOT}/.devcontainer/.deps-installed"

if [[ "${CODESPACES_RECOVERY_CONTAINER:-}" == "true" ]]; then
  echo "Recovery container detected. Skipping dependency installation."
  echo "Rebuild into the regular devcontainer to install project dependencies."
  exit 0
fi

if [[ -f "${MARKER_FILE}" ]]; then
  echo "Dependencies already installed. Skipping."
  python --version
  command -v java >/dev/null 2>&1 && java -version || echo "java not found"
  command -v nextflow >/dev/null 2>&1 && nextflow -version || echo "nextflow not found"
  exit 0
fi

run_privileged() {
  if command -v sudo >/dev/null 2>&1; then
    sudo "$@"
  else
    "$@"
  fi
}

retry() {
  local attempts="$1"
  shift

  local n=1
  while true; do
    if "$@"; then
      return 0
    fi

    if [[ "$n" -ge "$attempts" ]]; then
      echo "Command failed after ${attempts} attempts: $*" >&2
      return 1
    fi

    echo "Attempt ${n} failed. Retrying..." >&2
    n=$((n + 1))
  done
}

python --version

if ! command -v java >/dev/null 2>&1; then
  echo "Java not found. Installing OpenJDK 17 runtime..."
  if command -v apt-get >/dev/null 2>&1; then
    run_privileged apt-get update
    # Prefer Java 17+, but handle distro/package naming differences.
    JAVA_APT_PACKAGES=(
      openjdk-21-jre-headless
      openjdk-17-jre-headless
      default-jre-headless
    )

    java_pkg_installed="false"
    for pkg in "${JAVA_APT_PACKAGES[@]}"; do
      if run_privileged apt-get install -y --no-install-recommends "$pkg"; then
        java_pkg_installed="true"
        break
      fi
      echo "Package '$pkg' not available. Trying next Java package..."
    done

    if [[ "${java_pkg_installed}" != "true" ]]; then
      echo "Unable to install Java from apt-get using known package names." >&2
      exit 1
    fi
  elif command -v apk >/dev/null 2>&1; then
    run_privileged apk add --no-cache openjdk17-jre-headless
  else
    echo "No supported package manager found to install Java." >&2
    exit 1
  fi
fi
java -version

if ! command -v nextflow >/dev/null 2>&1; then
  echo "Nextflow not found. Installing Nextflow..."
  curl -s https://get.nextflow.io | bash
  run_privileged mv nextflow /usr/local/bin/nextflow
  run_privileged chmod +x /usr/local/bin/nextflow
fi
nextflow -version

if ! python -m pip --version >/dev/null 2>&1; then
  echo "pip not found for Python. Installing pip..."
  if command -v apt-get >/dev/null 2>&1; then
    run_privileged apt-get update
    run_privileged apt-get install -y --no-install-recommends python3-pip
  elif command -v apk >/dev/null 2>&1; then
    run_privileged apk add --no-cache py3-pip
  else
    echo "No supported package manager found to install pip." >&2
    exit 1
  fi
fi

PIP_BREAK_FLAG=()
if python -m pip help install 2>/dev/null | grep -q -- "--break-system-packages"; then
  PIP_BREAK_FLAG=(--break-system-packages)
fi

retry 3 python -m pip install "${PIP_BREAK_FLAG[@]}" --no-cache-dir --upgrade pip setuptools wheel

# Install torch and torchvision from the CPU wheels index to avoid GPU/CUDA resolver issues.
retry 3 python -m pip install "${PIP_BREAK_FLAG[@]}" --no-cache-dir --prefer-binary \
  --index-url https://download.pytorch.org/whl/cpu \
  "torch>=2.4,<2.7" \
  "torchvision>=0.19,<0.22"

LIBC_KIND="$(python - <<'PY'
import platform
print((platform.libc_ver()[0] or "").lower())
PY
)"

COMMON_PACKAGES=(
  "celltypist==1.7.1"
  "leidenalg==0.10.2"
  "scanpy==1.11.1"
  "scikit-learn==1.5.2"
)

retry 3 python -m pip install "${PIP_BREAK_FLAG[@]}" --no-cache-dir --prefer-binary \
  "${COMMON_PACKAGES[@]}"

if [[ "${LIBC_KIND}" == "musl" ]]; then
  echo "Skipping scvi-tools on musl/Alpine: no compatible jaxlib wheels available."
else
  retry 3 python -m pip install "${PIP_BREAK_FLAG[@]}" --no-cache-dir --prefer-binary \
    "scvi-tools==1.2.2"
fi

# Optional dependency: keep setup resilient if RAPIDS wheels are unavailable.
python -m pip install "${PIP_BREAK_FLAG[@]}" --no-cache-dir --prefer-binary rapids-singlecell || true

touch "${MARKER_FILE}"

echo "Post-create dependency install complete."

