FROM python:3.12

RUN pip install --no-cache torch torchvision --index-url https://download.pytorch.org/whl/cu126

RUN pip install --no-cache celltypist leidenalg scvi-tools scanpy scikit-learn rapids-singlecell
