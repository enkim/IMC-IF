# IMC-IF
Dual-modality imaging of immunofluorescence and imaging mass cytometry for high-resolution whole slide imaging with accurate single-cell segmentation


# Installation
conda create --name imc_process python=3.8 <br/> 
conda activate imc_process  <br/> 
pip install scimap <br/> 
python -m ipykernel install --user --name imc_process --display-name "imc_process" <br/> 
pip install scikit-image  <br/> 
pip install "dask[complete]" <br/> 
pip install 'napari[all]' <br/> 
pip install datashader <br/> 
conda install -c conda-forge holoviews <br/> 
pip install dask-image <br/> 
conda install -c conda-forge tifffile <br/> 
conda install -c conda-forge zarr <br/> 
pip install czifile <br/> 
