# Optional scVI/scANVI training parameters for nf-autoannotate.
# Values must be single settings, not parameter-sweep lists.

# Parameters for SCVI.setup_anndata.
param_input = {
    # "layer": "counts",
    # "batch_key": "sample_id",
    # "categorical_covariate_keys": ["donor", "chemistry"],
    # "continuous_covariate_keys": ["pct_counts_mt"],
}

# Parameters for scvi.model.SCVI.
model_input = {
    "n_layers": 2,
    "n_latent": 20,
    # "n_hidden": 128,
    # "dropout_rate": 0.1,
}

# Parameters for SCVI.train.
train_input = {
    "max_epochs": 100,
    "batch_size": 512,
    "early_stopping": True,
    "early_stopping_patience": 10,
}

# Optional parameters for SCANVI.train. If omitted, train_input is reused.
scanvi_train_input = {
    "max_epochs": 100,
    "batch_size": 512,
    "early_stopping": True,
    "early_stopping_patience": 10,
}
