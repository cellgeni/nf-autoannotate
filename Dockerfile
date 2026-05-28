FROM python:3.12

RUN pip install --no-cache-dir torch torchvision --index-url https://download.pytorch.org/whl/cu126

RUN pip install --no-cache-dir celltypist leidenalg scvi-tools scanpy scikit-learn==1.7.2

RUN pip install --no-cache-dir "cudf-cu12==26.2.*"
