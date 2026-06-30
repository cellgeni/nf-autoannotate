# nf-scVI-metrics-style input file for nf-autoannotate compatibility mode.
# nf-autoannotate accepts exactly one value per model/training parameter.
# Single-item lists from nf-scVI-metrics are allowed and unwrapped; multi-value
# sweep lists are rejected.

# Dataset selectors from nf-scVI-metrics, such as adata_path, are ignored here
# because nf-autoannotate uses --ref_h5ad and --query_h5ad.
model_input = {
    "n_latent": [30],
    "n_layers": [2],
    "dispersion": ["gene"],
    "gene_likelihood": ["nb"],
    "max_epochs": 100,
    "batch_size": 512,
    "early_stopping": True,
    "early_stopping_patience": 10,
    "counts_layer": None,
}

# Parameters for SCVI.setup_anndata.
param_input = {
    "batch_key": "sample_id",
    "continuous_covariate_keys": [],
    "categorical_covariate_keys": [],
}

# Optional nf-autoannotate extension for the scANVI fine-tuning step.
# If omitted, model_input training parameters are reused.
# scanvi_train_input = {
#     "max_epochs": 20,
#     "batch_size": 512,
# }

# Optional nf-autoannotate extension for query-mapping SCANVI.train.
# query_train_input = {
#     "max_epochs": 100,
#     "batch_size": 512,
#     "check_val_every_n_epoch": 10,
#     "plan_weight_decay": 0.0,
# }
