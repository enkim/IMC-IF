# IMC-IF
Dual-modality imaging of immunofluorescence and imaging mass cytometry for high-resolution whole slide imaging with accurate single-cell segmentation


# Installation
conda create --name imc_process python=3.8
conda activate imc_process
pip install scimap
python -m ipykernel install --user --name imc_process --display-name "imc_process"
pip install scikit-image 
pip install "dask[complete]"
pip install 'napari[all]'
pip install datashader
conda install -c conda-forge holoviews
pip install dask-image 
conda install -c conda-forge tifffile
conda install -c conda-forge zarr 
pip install czifile 
