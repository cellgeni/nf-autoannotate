"""Check that main.nf only passes flags the Python scripts actually accept.

The stubbed pipeline test replaces every `script:` block with a fake, so the real
command lines are never executed by any test. A stale flag left behind after a script's
argparse changed would therefore pass the whole suite and only fail in a real container
run, minutes into a job.

This closes that gap statically: it reads the flags out of each process in main.nf and
compares them against the arguments each script declares. Both sides are parsed as text
and AST, so nothing here needs scvi-tools, celltypist, or a container.
"""

import ast
import re
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
MAIN_NF = REPO / "main.nf"
SCRIPTS = REPO / "scripts"

PROCESS_PATTERN = re.compile(r"^process (\w+) \{$", re.MULTILINE)
SCRIPT_INVOCATION_PATTERN = re.compile(r"\$\{projectDir\}/scripts/(\w+\.py)")
FLAG_PATTERN = re.compile(r"(?<![\w-])--([a-z0-9][a-z0-9-]*)")


def process_blocks() -> dict[str, str]:
    """Split main.nf into {process name: body} using the top-level process declarations."""
    source = MAIN_NF.read_text(encoding="utf-8")
    matches = list(PROCESS_PATTERN.finditer(source))
    blocks = {}
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(source)
        blocks[match.group(1)] = source[match.start():end]
    return blocks


def script_command(body: str, script_name: str) -> str:
    """Return just the shell command, from the python3 invocation to its closing quotes.

    Bounded deliberately: the final process block runs to the end of the file, so an
    unbounded slice would sweep up error messages in the workflow body that mention
    parameters like --tools.
    """
    start = body.index(f"scripts/{script_name}")
    end = body.index('"""', start)
    return body[start:end]


def declared_arguments(script_name: str) -> set[str]:
    """Return the long-form flags a script's argparse accepts."""
    tree = ast.parse((SCRIPTS / script_name).read_text(encoding="utf-8"))
    flags = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not isinstance(node.func, ast.Attribute) or node.func.attr != "add_argument":
            continue
        for arg in node.args:
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str) and arg.value.startswith("--"):
                flags.add(arg.value[2:])
    return flags


class CliContractTests(unittest.TestCase):
    def test_every_process_only_passes_flags_its_script_accepts(self):
        checked = 0
        for process_name, body in process_blocks().items():
            scripts = set(SCRIPT_INVOCATION_PATTERN.findall(body))
            if not scripts:
                continue
            self.assertEqual(
                len(scripts), 1, f"{process_name} invokes more than one script: {sorted(scripts)}"
            )
            script_name = scripts.pop()

            # Flags appear literally in the command even inside Groovy conditionals.
            passed = set(FLAG_PATTERN.findall(script_command(body, script_name)))
            accepted = declared_arguments(script_name)

            unknown = sorted(passed - accepted)
            self.assertEqual(
                unknown,
                [],
                f"{process_name} passes flags {script_name} does not accept: {unknown}. "
                f"{script_name} accepts: {sorted(accepted)}",
            )
            checked += 1

        # Guards against the regexes silently matching nothing and the test passing vacuously.
        self.assertGreaterEqual(checked, 6, "expected to check at least six script-invoking processes")

    def test_required_script_arguments_are_always_supplied(self):
        """A required argument must be passed unconditionally, not inside a Groovy ternary."""
        for process_name, body in process_blocks().items():
            scripts = set(SCRIPT_INVOCATION_PATTERN.findall(body))
            if not scripts:
                continue
            script_name = scripts.pop()
            command = script_command(body, script_name)

            tree = ast.parse((SCRIPTS / script_name).read_text(encoding="utf-8"))
            required = set()
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                if not isinstance(node.func, ast.Attribute) or node.func.attr != "add_argument":
                    continue
                is_required = any(
                    keyword.arg == "required"
                    and isinstance(keyword.value, ast.Constant)
                    and keyword.value.value is True
                    for keyword in node.keywords
                )
                if not is_required:
                    continue
                for arg in node.args:
                    if isinstance(arg, ast.Constant) and isinstance(arg.value, str) and arg.value.startswith("--"):
                        required.add(arg.value[2:])

            for flag in sorted(required):
                with self.subTest(process=process_name, flag=flag):
                    self.assertRegex(
                        command,
                        rf"(?<![\w-])--{re.escape(flag)}(?![\w-])",
                        f"{process_name} never passes required flag --{flag} of {script_name}",
                    )


if __name__ == "__main__":
    unittest.main()
