# conda activate deepcell-env 
# IMC_IF_11 36marker_IFtile_mesmer_DAPI_Ecadherin 

from deepcell.applications import Mesmer
import numpy as np
import os
from PIL import Image 
import pandas as pd
from skimage.util import img_as_uint
from skimage.io import imread,imsave,imshow
from skimage.color import label2rgb

os.chdir('/home/groups/ChangLab/kimeunn/CRUK/new_July_2021_EN') # CRC project box folder 

 # Create the application
app = Mesmer()

df=pd.read_csv('sample/IFtilelist.csv')

for index, row in df.iterrows():
    if row['dapimax'] > 10 :
        imgpath='IFtile/'+ row['tile'].split('/')[-1].replace('_mask.png','')
        filename=imgpath.split('/')[-1].replace('_y_','_y_0').replace('y_010','y_10').replace('y_011','y_11')+"_mask.png"
        
        # Load the images
        im1 = imread(imgpath)[:,:,0] # input DAPI channel  
        im2 = imread(imgpath)[:,:,2] # input Ecad channel  

        # Combined together and expand to 4D 
        im = np.stack((im1, im2), axis=-1)
        im = np.expand_dims(im,0)
        
        labeled_image = app.predict(im, image_mpp=0.325) # input your image's resolution, ex. CRC: 0.325um/pix  
        mask=labeled_image[0,:,:,0] 
        Image.fromarray(labeled_image).save(f'IFtilemask_mesmer/{filename}')
        print(filename)
    else:
        imgpath='IFtile/'+ row['tile'].split('/')[-1].replace('_mask.png','')
        filename=imgpath.split('/')[-1].replace('_y_','_y_0').replace('y_010','y_10').replace('y_011','y_11')+"_mask.png"
        or_img=imread(imgpath)
        empty_img=np.zeros([or_img.shape[0], or_img.shape[1]],dtype=np.uint16)
        Image.fromarray(empty_img).save(f'IFtilemask_mesmer/{filename}')
        print(filename)
        
print('imgmask done')
