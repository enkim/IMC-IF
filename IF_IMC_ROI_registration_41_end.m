% IF - IMC ROI registration 
% New dataset
clear all;
clc;
close all;

% image path
HE_path = './1. Brightfield/3. HE/';
IF_path = './2. IF/';
IMC_path = './3. TIFF/';
output_path = './registered_IMC/';
mkdir(output_path);



HE_filename = dir(sprintf('%s*.tif', HE_path));
IF_filename = dir(sprintf('%s*.tif', IF_path));


% parameter
HE_scan = 264.5833*10^-3; % umm/pxl 
IMC_scan = 1.0; % umm/pxl
IF_scan = 0.325; % umm/pxl
N_smpl = 2000; % num. of feature pts.
    
%% image reading
%I_HE = imread(sprintf('%s%s',HE_path,HE_filename.name));
I_IF = imread(sprintf('%s%s',IF_path, IF_filename.name));


% image resize
I_IF_resize = imresize(I_IF, IF_scan/IMC_scan);
I_IF_resize = imrotate(I_IF_resize,180);

%I_HE_resize = imresize(I_HE, HE_scan/IMC_scan); 
%I_HE_resize = imrotate(I_HE_resize,180);

% reference feature
IF_DAPI = I_IF_resize(:,:,1);
IF_CK = I_IF_resize(:,:,2);
IF_CD45 = I_IF_resize(:,:,3);


I_ref = imbinarize(IF_DAPI, graythresh(IF_DAPI)); % reference IF DAPI
ptsRef_all = detectSURFFeatures(I_ref);
%ptsRef=ptsRef_all.selectStrongest(N_smpl*5);


%%
% ROI
roi_name = dir(sprintf('%sROI*', IMC_path));

marker_name = dir(sprintf('%s%s/*.tiff', IMC_path,roi_name(1).name));

for ch=41:length(marker_name)
    %ch = 44;
    
    fprintf('marker: %s\n', marker_name(ch).name);
    
    
    Registered_IMC = zeros( size(I_ref));
    for r=1:length(roi_name)
      

        fprintf('ROI: %d\n', r);

        % find DNA image
        I_IMC_DNA = [];
        IMC_DNA_name = dir(sprintf('%s%s/*DNA.ome.tiff', IMC_path, roi_name(r).name));
        I_IMC_DNA(:,:,1) = imread(sprintf('%s%s/%s', IMC_path, roi_name(r).name, IMC_DNA_name(1).name));
        I_IMC_DNA(:,:,2) = imread(sprintf('%s%s/%s', IMC_path, roi_name(r).name, IMC_DNA_name(2).name));
        I_obj = max(I_IMC_DNA(:,:,1), I_IMC_DNA(:,:,2));

        I_obj = imadjust(uint8(I_obj/max(I_obj(:))*255));
        I_obj = imbinarize(I_obj,graythresh(I_obj));

        
        % read IMC image
        I_IMC = [];
        I_IMC = imread(sprintf('%s%s/%s', IMC_path, roi_name(r).name,marker_name(ch).name));

        % correlation based
        c = normxcorr2( I_obj, I_ref);
        %figure, surf(c), shading flat
        [max_c, imax] = max(abs(c(:)));
        [ypeak, xpeak] = ind2sub(size(c),imax(1));
        corr_offset = [(xpeak-size(I_obj,2)) 
                   (ypeak-size(I_obj,1))];

        % total offset
        offset = corr_offset;
        xoffset = offset(1);
        yoffset = offset(2);

        xbegin = round(xoffset+1);
        xend   = round(xoffset+ size(I_obj,2));
        ybegin = round(yoffset+1);
        yend   = round(yoffset+size(I_obj,1));

        registered_obj = (zeros(size(I_ref)));
        registered_obj(ybegin:yend, xbegin:xend,:) = I_obj;
        registered_obj = registered_obj > 0;

        registered_IMC = (zeros(size(I_ref)));
        registered_IMC(ybegin:yend, xbegin:xend,:) = I_IMC;
       % imagesc(imfuse(I_ref, registered_obj, 'falsecolor', 'Scaling','joint', 'ColorChannels', [1 2 0]));


        %

        ptsObj_all = []; ptsObj = [];
        ptsObj_all = detectSURFFeatures(registered_obj);
        %ptsObj =ptsObj_all.selectStrongest(N_smpl);


        xv = [ xbegin xbegin xend xend]';
        yv = [ yend ybegin ybegin yend]';
        % select feature within location

        ptsObj = ptsObj_all( inpolygon( ptsObj_all.Location(:,1), ptsObj_all.Location(:,2), xv, yv));
        ptsRef = ptsRef_all( inpolygon( ptsRef_all.Location(:,1), ptsRef_all.Location(:,2), xv, yv));



        [featuresRef, validPtsRef] = extractFeatures( I_ref, ptsRef);
        [featuresObj, validPtsObj] = extractFeatures( registered_obj, ptsObj);

        indxPairs = matchFeatures(featuresRef, featuresObj);
        matchedRef = validPtsRef(indxPairs(:,1));
        matchedObj = validPtsObj(indxPairs(:,2));

        [tform, inlierDistorted, inlierOriginal,status] = estimateGeometricTransform(...
              matchedObj, matchedRef, 'similarity', 'MaxNumTrials',500, 'Confidence',99.00, 'MaxDistance',5);

%          I_over = imfuse(I_ref, imwarp(registered_obj, tform, 'OutputView', imref2d(size(I_ref))), ...
%                              'falsecolor', 'Scaling','joint', 'ColorChannels', [1 2 0]);


        registered_IMC = imwarp(registered_IMC, tform, 'OutputView', imref2d(size(I_ref)));
        
        
        % Using max projection
        Registered_IMC = max( registered_IMC, Registered_IMC);
                         
        %imagesc(I_over);
        fprintf('ROI %d registration done\n', r);
    end
    
    imwrite( imrotate(uint16(Registered_IMC),-180), sprintf('%sRegistered_%s', output_path, marker_name(ch).name), 'Tiff');
    
    % checking
    % imagesc(imfuse(I_ref, imbinarize(Registered_IMC, graythresh(Registered_IMC)), 'falsecolor', 'Scaling','joint', 'ColorChannels', [1 2 0]));
end
