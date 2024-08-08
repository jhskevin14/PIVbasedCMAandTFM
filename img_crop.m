function [t_num, crpdimg_dir]=img_crop(t0, t_num, img_dir, domain, im_ph_fst)
%% ===================== Center of mass in the domain =====================
s = regionprops(domain,'centroid', 'Area');

[~,index] = max([s.Area]);
cent = cat(1, s(index).Centroid);
cent_r = round(cent);

%%  ============================= Finding ROI =============================

[M, ~] = size(im_ph_fst);
wth = (floor(M/32) - 2)*32;

i1 = cent_r(2) - wth/2+1;  
i2 = cent_r(2) + wth/2;

if i1<=10
    md = round((M-i2)/2)+i1;
    i2 = i2-i1+md;
    i1 = md;
end

j1 = cent_r(1) - wth/2+1;
j2 = cent_r(1) + wth/2;

%%  ================= Dedrifting and cropping all images =================
disp(['======== Cropping images (' num2str(t0) 'of ' num2str(t_num) ') ========']);
im_bd_fst = imread([img_dir '/Rawdata/bead.tif'],t0);

im_ph_fst = im_ph_fst(i1:i2, j1:j2);
im_bd_fst = im_bd_fst(i1:i2, j1:j2);

mkdir([img_dir '/Cropped_img']);
crpdimg_dir = [img_dir, '/Cropped_img'];

if exist([crpdimg_dir '/phase_crpt.tif'],'file') == 2
    delete([crpdimg_dir '/phase_crpt.tif']);
end
if exist([crpdimg_dir '/bead_crpt.tif'],'file') == 2
    delete([crpdimg_dir '/bead_crpt.tif']);
end

imwrite(im_ph_fst,[crpdimg_dir '/phase_crpt.tif'],'tif','writemode','append');
imwrite(im_bd_fst,[crpdimg_dir '/bead_crpt.tif'],'tif','writemode','append');

for t = t0+1:t_num
    disp(['======== Cropping images (' num2str(t) ' of ' num2str(t_num) ') ========']);
    im_bd_cur = imread([img_dir '/Rawdata/bead.tif'],t);
    im_bd_tmp = im_bd_cur(i1:i2, j1:j2);
    [~, ~, shx, shy, ~] = disp_on_blocks_fast(double(im_bd_fst), ...   % reference image
                          double(im_bd_tmp), ...                      % current image
                          size(im_bd_fst,1), ...                       % image size
                          size(im_bd_fst,2), ...                       % image size
                          0, ...                                       % overlap between blocks
                          0, ...                                       % padding for the PIV, usually 0
                          'hanning', ...                               % window
                          '2DCentroid', ...                            % Subpixel interpolation method
                          0, ...                                       % threshold of center of mass calculation
                          10, ...                                      % Size of the window
                          4);                                          % number of iterations
    shx_int = round(shx);
    shy_int = round(shy);
    
    im_bd_cur = im_bd_cur(i1+shy_int:i2+shy_int, j1+shx_int:j2+shx_int);
    
    im_ph_cur = imread([img_dir '/Rawdata/phase.tif'],t);
    im_ph_cur = im_ph_cur(i1+shy_int:i2+shy_int, j1+shx_int:j2+shx_int);
    
    imwrite(im_ph_cur,[crpdimg_dir '/phase_crpt.tif'],'tif','writemode','append');
    imwrite(im_bd_cur,[crpdimg_dir '/bead_crpt.tif'],'tif','writemode','append');
end
disp(['======== Cropping trypsin image ========']);
im_bd_tr = imread([img_dir '/Rawdata/bead_tr.tif']);
im_bd_tr2 = im_bd_tr(i1:i2, j1:j2);
[~, ~, shx, shy, ~] = disp_on_blocks_fast(double(im_bd_fst), ...   % reference image
                      double(im_bd_tr2), ...                        % current image
                      size(im_bd_fst,1), ...                       % image size
                      size(im_bd_fst,2), ...                       % image size
                      0, ...                                       % overlap between blocks
                      0, ...                                       % padding for the PIV, usually 0
                      'hanning', ...                               % window
                      '2DCentroid', ...                            % Subpixel interpolation method
                      0, ...                                       % threshold of center of mass calculation
                      10, ...                                      % Size of the window
                      4);                                          % number of iterations
 shx = round(shx);
 shy = round(shy);
 
 im_bd_tr = im_bd_tr(i1+shy:i2+shy, j1+shx:j2+shx);
 imwrite(im_bd_tr,[crpdimg_dir '/bead_tr_crpt.tif'],'tif');
 
end
 