#!/usr/bin/env python3

"""Rewrite main.nf with fake annotation steps so the channel plumbing can be tested.

The real CellTypist, scANVI, and PanHumanPy steps need GPUs and containers, so a test
that exercises the workflow end to end has to replace them. Everything else - sample-ID
resolution, per-query fan-out, grouping, joining, and publishing - is the real pipeline.

Usage: annotation_stubs.py <main.nf> <output main.nf>
"""

import sys
from pathlib import Path


STUBS = {
    "PANHUMANPY_ANNOTATE": '''    script:
    """
    python3 -c "
import anndata as ad, pandas as pd
q = ad.read_h5ad('${query_h5ad}')
pd.DataFrame({'cell_id': q.obs_names.astype(str),
 'panhumanpy_full_hierarchical_label': 'A|B', 'panhumanpy_level_zero_label': 'A',
 'panhumanpy_final_level_label': 'B', 'panhumanpy_confidence': 0.9}).to_csv('panhumanpy_predictions.csv', index=False)
"
    """
}
''',
    "CELLTYPIST_TRAIN_MODEL": '''    script:
    """
    echo stub > celltypist_model.pkl
    echo '{}' > celltypist_model_metadata.json
    """
}
''',
    "CELLTYPIST_ANNOTATE": '''    script:
    """
    python3 -c "
import anndata as ad, pandas as pd
q = ad.read_h5ad('${query_h5ad}')
pd.DataFrame({'cell_id': q.obs_names.astype(str),
 'celltypist_predicted_label': 'B', 'celltypist_majority_voting': 'B',
 'celltypist_confidence': 0.8}).to_csv('celltypist_predictions.csv', index=False)
"
    """
}
''',
    "SCANVI_TRAIN_MODEL": '''    script:
    """
    mkdir -p scvi_model scanvi_model
    echo '{}' > scvi_model_metadata.json
    echo '{}' > scanvi_model_metadata.json
    """
}
''',
    "SCANVI_ANNOTATE": '''    script:
    """
    python3 -c "
import anndata as ad, numpy as np, pandas as pd
q = ad.read_h5ad('${query_h5ad}')
pd.DataFrame({'cell_id': q.obs_names.astype(str),
 'scanvi_label': 'B', 'scanvi_confidence': 0.7}).to_csv('scanvi_predictions.csv', index=False)
latent = ad.AnnData(obs=pd.DataFrame(index=q.obs_names.astype(str)))
latent.obsm['X_scanvi'] = np.zeros((q.n_obs, 2))
latent.write_h5ad('scanvi_latent.h5ad')
"
    """
}
''',
}


def stub_main_nf(source: str) -> str:
    for name, stub in STUBS.items():
        start = source.index(f"process {name} {{")
        script_start = source.index("    script:", start)
        end = source.index("\n}\n", script_start) + len("\n}\n")
        source = source[:script_start] + stub + source[end:]
    return source


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    source = Path(sys.argv[1]).read_text(encoding="utf-8")
    Path(sys.argv[2]).write_text(stub_main_nf(source), encoding="utf-8")


if __name__ == "__main__":
    main()
