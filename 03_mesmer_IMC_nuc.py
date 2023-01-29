# conda activate deepcell-env 

from skimage.io import imread
from deepcell.applications import Mesmer
import numpy as np
import os
import glob
from PIL import Image 
from skimage.util import img_as_uint
from skimage.io import imread,imsave,imshow
from skimage import img_as_ubyte

os.chdir('/home/groups/ChangLab/kimeunn/CRUK/new_July_2021_EN') # CRC project box folder 

# Create the application
app = Mesmer()

samplenames=sorted([filename.split('/')[-1].split('.')[0].split('_')[-1] for filename in glob.glob('IMC_registration/*')])
n=0 

for n in range(len(samplenames)): 
    samplename= samplenames[n]
    DNA_path=glob.glob(f'IMC_registration/*{samplename}*/*191_DNA*')[0]
    Ecad_path=glob.glob(f'IMC_registration/*{samplename}*/*E-Ca*')[0]

    # Load the images
    im1 = np.array(imread(DNA_path),np.uint32) # input DNA IMC  
    im2 = np.array(imread(DNA_path),np.uint32) # input DNA IMC  

    # Combined together and expand to 4D 
    im = np.stack((im1, im2), axis=-1)
    im = np.expand_dims(im,0)

    # create the lab
    labeled_image = app.predict(im, image_mpp=1.0) # input your image's resolution, ex. CRC: 0.325um/pix  
    mask=labeled_image[0,:,:,0] 
    mask=np.array(mask,np.uint32)
    Image.fromarray(mask).save(f'IMC_mesmer_nuc/{samplename}.tiff')


    print('cell number after resizing:', np.max(mask))
    print('image type after resizing:',mask.dtype)

print('all done')