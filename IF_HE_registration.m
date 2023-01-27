% IF - H&E ROI registration 
% New dataset
clear all;
clc;
close all;

% image path
HE_path = './1. Brightfield/3. HE/';

% temporarily I used H&E path

%%%%% need to update %%%%%%%%%%%%
% insert : mask image path in HE_path   
HE_mask_path = HE_path; % H&E mask image directory (* need to update) 
%%%%% need to update %%%%%%%%%%%%


IF_path = './2. IF/';
IMC_path = './3. TIFF/';
output_path = './registered_HE/';
mkdir(output_path);


HE_filename = dir(sprintf('%s*.tif', HE_path));
HE_mask_filename = dir(sprintf('%s*.tif', HE_mask_path)); 
IF_filename = dir(sprintf('%s*.tif', IF_path));


% parameter
HE_scan = 0.22; % umm/pxl 
IMC_scan = 1.0; % umm/pxl
IF_scan = 0.325; % umm/pxl



%% image reading
I_HE = imread(sprintf('%s%s',HE_path,HE_filename.name));
I_HE_mask = imread(sprintf('%s%s',HE_mask_path,HE_mask_filename.name));
I_IF = imread(sprintf('%s%s',IF_path, IF_filename.name));


% image resize
I_IF_resize = imresize(I_IF, IF_scan/IMC_scan);
I_IF_resize = imrotate(I_IF_resize,180);

I_HE_resize = imresize(I_HE, HE_scan/IMC_scan); %HE_scan/IF_scan); 
I_HE_resize = imrotate(I_HE_resize,180);

I_HE_mask_resize = imresize(I_HE_mask, HE_scan/IMC_scan);
I_HE_mask_resize = imrotate(I_HE_mask_resize,180);



%% registration

% reference image
IF_DAPI = I_IF_resize(:,:,1);
I_ref = imbinarize(IF_DAPI, graythresh(IF_DAPI)*0.5); % reference IF DAPI
%ptsRef_all = detectSURFFeatures(I_ref);
%ptsRef=ptsRef_all.selectStrongest(N_smpl);


I_HE_gray = (rgb2gray(I_HE_resize));
I_HE_obj = imcomplement(imbinarize( I_HE_gray, graythresh(I_HE_gray)*0.65));

I_HE_BD = imfill(imerode(I_HE_gray > 0, strel('disk',20)), 'holes');

I_HE_obj( find(I_HE_BD == 0)) = 0;

I_obj = I_HE_obj;

%%
%%
% correlation based

%I_ref_org = I_ref;
%I_ref = zeros(floor(size(I_obj,1)*1.25), floor(size(I_obj,2)*1.25))>0; % make sure reference is bigger than obj
%I_ref(1:size(I_ref_org,1), 1:size(I_ref_org,2)) = I_ref_org;


c = normxcorr2( I_obj, I_ref);
%figure, surf(c), shading flat
[max_c, imax] = max(abs(c(:)));
[ypeak, xpeak] = ind2sub(size(c),imax(1));
corr_offset = [(xpeak-size(I_obj,2)) 
           (ypeak-size(I_obj,1))];

% total offset
offset = corr_offset;
xoffset = max(offset(1),0);
yoffset = max(offset(2),0);

xbegin = round(xoffset+1);
xend   = round(xoffset+ size(I_obj,2));
ybegin = round(yoffset+1);
yend   = round(yoffset+size(I_obj,1));

registered_obj = (zeros(size(I_ref)));
registered_obj(ybegin:yend, xbegin:xend,:) = I_obj;
registered_obj = registered_obj > 0;

registered_mask = (zeros(size(I_ref)));
registered_mask(ybegin:yend, xbegin:xend,:) = I_HE_mask_resize;


registered_HE = zeros(size(I_ref,1),size(I_ref,2),3);
for ch=1:3
    registered_HE(ybegin:yend, xbegin:xend,ch)= I_HE_resize(:,:,ch);
end


I_over = imfuse(I_ref, registered_obj, ...
                              'falsecolor', 'Scaling','joint', 'ColorChannels', [1 2 0]);
figure; imagesc(I_over)


%%
scale = 0.25;
I_ref_scale = imresize(I_ref,scale); % temp
registered_obj_scale = imresize(registered_obj, scale); % temp


moving_scale = uint8(registered_obj_scale);
fixed_scale = uint8(I_ref_scale);
[optimizer,metric] = imregconfig('multimodal');
%movingRegisteredDefault = imregister(moving,fixed,'similar',optimizer,metric);

optimizer.MaximumIterations = 500;
optimizer.InitialRadius = 3.5e-3;
%optimizer.InitialRadius = optimizer.InitialRadius/3;


%[movingRegisteredAdjustedInitialRadius, tform] = imregister(moving,fixed,'rigid',optimizer,metric);

tform_scale = imregtform(moving_scale, fixed_scale, 'similar', optimizer, metric);

movingReg_scale = imwarp(moving_scale, tform_scale, 'OutputView',imref2d(size(fixed_scale)));


tform = tform_scale;

tform.T(3,1) = tform_scale.T(3,1)/scale;
tform.T(3,2) = tform_scale.T(3,2)/scale;


fixed = uint8(I_ref);
moving = uint8(registered_obj);
movingReg = imwarp(moving, tform, 'OutputView', imref2d(size(fixed)));

I_over2 = imfuse( fixed, movingReg, ...
                              'falsecolor', 'Scaling','joint', 'ColorChannels', [1 2 0]);

figure; imagesc(I_over2)

reg_mask = imwarp(registered_mask, tform, 'OutputView', imref2d(size(fixed)));

registered_HE_final = zeros(size(I_ref,1),size(I_ref,2),3);
for ch=1:3
    reg_HE_final(:,:,ch) = imwarp(registered_HE(:,:,ch), tform, 'OutputView', imref2d(size(fixed)));
end


%imwrite( imrotate(uint8(reg_mask),-180), sprintf('%sRegistered_HE_mask_%s', output_path,HE_filename.name), 'Tiff');

imwrite( imrotate(uint8(reg_HE_final),-180), sprintf('%sRegistered_HE_%s', output_path,HE_filename.name), 'Tiff');
    

