%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMC / HE / IF registration 
% Chang Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
clc;
close all;

smpl_name = dir('./ROI*');

output_dir = './output/';
%mkdir(output_dir);

mask_dir = './IFmask/'; % input
mask_out_dir = './reg_mask/'; % output
%mkdir(mask_out_dir);

CD45IF_out_dir= './reg_CD45IF/';
%mkdir(CD45IF_out_dir);


for smpl=1:length(smpl_name)
    close all

file loading
file folder / name
    fname = smpl_name(smpl).name;
    in_dir = sprintf('./%s/', fname);
    fprintf('%s', in_dir);
    
    IMC_DNA_fname = dir(sprintf('%s*DNA.ome.tiff', in_dir));    

    HE_name = dir(sprintf('%s*_HE.tif', in_dir));
    IF_DAPI_name = dir(sprintf('%s*_c1_*.tif', in_dir));
    IF_CD45_name = dir(sprintf('%s*_c2_*.tif', in_dir));

    HE_scan = 264.5833*10^-3; % umm/pxl 
    IMC_scan = 1.0; % umm/pxl
    IF_scan = 0.325; % umm/pxl
    N_smpl = 1000;


    % read IMC data
    I_IMC_DNA = uint16([]);
    for i=1:length(IMC_DNA_fname)
        I_IMC_DNA(:,:,i) = imread(sprintf('%s%s', in_dir, IMC_DNA_fname(i).name));
    end

    I_IMC_DNA_max = max(I_IMC_DNA(:,:,1), I_IMC_DNA(:,:,2)); % max projection
    I_IMC_DNA_max = uint8( imadjust(I_IMC_DNA_max)/255);
    %I_IMC_DNA_max = imgaussfilt(I_IMC_DNA_max,0.5);



    % read H&E
%     I_HE = imread(sprintf('%s%s', in_dir, HE_name.name));
%     I_HE_res =  imresize(I_HE, HE_scan/IMC_scan);
%     I_HE_gray = imcomplement(rgb2gray(I_HE_res));



    % read IF data
    IF_DAPI = imread(sprintf('%s%s', in_dir, IF_DAPI_name.name));
    IF_DAPI_res = imresize(IF_DAPI, IF_scan/IMC_scan);
    IF_CD45 = imread(sprintf('%s%s', in_dir, IF_CD45_name.name));
    IF_CD45_res = imresize(IF_CD45, IF_scan/IMC_scan);
    
    % read segment mask
    IF_mask = imread(sprintf('%s%s_mask.png', mask_dir, fname)); 
    IF_mask_res = imresize(IF_mask, IF_scan/IMC_scan,'nearest');

Registration
    I_ref = imbinarize( I_IMC_DNA_max, graythresh(I_IMC_DNA_max));
    ptsRef_all =detectSURFFeatures(I_ref);
    ptsRef=ptsRef_all.selectStrongest(N_smpl);

    for task=1:1
        
        task = 2; % IF -> IMC only
        switch(task)
            case 1,
                I_obj = imbinarize( I_HE_gray, graythresh(I_HE_gray)); % HE -> IMC
                scale_resize = HE_scan/IMC_scan;
            case 2, 
                I_obj = imbinarize( IF_DAPI_res, graythresh(IF_DAPI_res)); % IF -> IMC
                scale_resize = IF_scan/IMC_scan;
        end

        
        ptsObj_all = []; ptsObj = [];
        ptsObj_all = detectSURFFeatures(I_obj);
        ptsObj =ptsObj_all.selectStrongest(N_smpl);


        [featuresRef, validPtsRef] = extractFeatures( I_ref, ptsRef);
        [featuresObj, validPtsObj] = extractFeatures( I_obj, ptsObj);

        indxPairs = matchFeatures(featuresRef, featuresObj);
        matchedRef = validPtsRef(indxPairs(:,1));
        matchedObj = validPtsObj(indxPairs(:,2));

       [tform_scale, inlierDistorted, inlierOriginal,status] = estimateGeometricTransform(...
          matchedObj, matchedRef, 'similarity', 'MaxNumTrials',5000, 'Confidence',99.99, 'MaxDistance',2);

        tform = [];
        tform = tform_scale;
        tform.T(3,1) = tform_scale.T(3,1)/scale_resize;
        tform.T(3,2) = tform_scale.T(3,2)/scale_resize;

       % outputView_scale = imref2d(size(I_ref)); % outputView for scaled one


        I_over_wo = imfuse(I_ref, I_obj, ...
                         'falsecolor', 'Scaling','joint', 'ColorChannels', [1 2 0]);

        I_over = imfuse(I_ref, imwarp(I_obj, tform_scale, 'OutputView', imref2d(size(I_ref))), ...
                         'falsecolor', 'Scaling','joint', 'ColorChannels', [1 2 0]);

        IF_CD45_reg = imwarp(IF_CD45_res, tform_scale, 'nearest', 'OutputView', imref2d(size(I_ref)));
        
        imwrite(IF_CD45_reg, sprintf('%s%s_reg_CD45.png', CD45IF_out_dir, fname), 'png');
          
        
        figure('pos',[10 10 1600 800]);
        subplot(121); imagesc(I_over_wo); title('before registration');
        subplot(122); imagesc(I_over); title('after registration');

        saveas(gcf,sprintf('%sregistered_%s.png', output_dir,fname));
        saveas(gcf,'filename.png')

        
        
    end

end