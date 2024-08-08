clear;
close all;
clc;

% ========================= 0_Basic input setup ===========================
th = 100;         % Gel thickness
E = 6e3;          % Young's modulus
nu = 0.48;        % Poisson's ratio
pix_size = 1.14; % 1.14;  % um (JuLi with 4x) (0.88 pixels/um)
int_time = 10; % Interval of timelaps (min)

grid_resolution_ph = 64;    % Resolution of the PIV grid for phase images
grid_space_ph = 16;
overlap_ph = 1-grid_space_ph/grid_resolution_ph;  % Overlapping ratio between subset blocks for phase images

grid_resolution_bd = 64;    % Resolution of the PIV grid for bead images
grid_space_bd = 16;
overlap_bd = 1-grid_space_bd/grid_resolution_bd;  % Overlapping ratio between subset blocks for bead images
% =========================================================================

% ============== 1_Addressing directories (codes and data) ===============

% disp('======== Select the folder containing these Matlab codes ========');
% code_dir = uigetdir('','Select the folder containing these Matlab codes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% manually %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
code_dir = 'C:\Users\jhske\Documents\Projects\PIVbasedCMAandTFM';
% code_dir = '/Users/jang_hwanseok/Dropbox/MatlabCode_HS/161201_codeHS';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(code_dir);

%% main directory change
load('inp_setup.mat');

tf = strcmp(pwd,img_dir);

if tf == 0
    img_dir = pwd;
    crpdimg_dir = [img_dir '/Cropped_img'];
end
code_dir = 'C:\Users\jhske\Documents\Projects\PIVbasedCMAandTFM';
% code_dir = '/Users/jang_hwanseok/Dropbox/MatlabCode_HS/161201_codeHS';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(code_dir);

save('inp_setup.mat', 'img_dir', 'crpdimg_dir', 'code_dir','-append');

%% ==================== 2_JuLI raw data to multi tiff =====================

% Samplestyle = menu('Sample style', 'Samples into individual positions', 'Samples needed stitching positions');
% switch Samplestyle
%     case 1
        [img_dir]=juliraw_to_multitif;
%     case 2
%         [img_dir]=juliraw_to_stitched_multitif;
% end 
cd(img_dir);
%% =========================== 3_Cropping images ==========================
t = 1;
t_num = length(imfinfo([img_dir '/Rawdata/phase.tif']));

im_ph_fst = imread([img_dir '/Rawdata/phase.tif'],t);

filter_th = 0.14; % Threshold for Matlab Sobel filter
dilate_rad = 4; % Radius of circular kernal used to dilate bright spots

modif = 0;
[M, N] = size(im_ph_fst);
[~, xedge, yedge, domain] = find_boundary(im_ph_fst, filter_th, dilate_rad);

sldr_1.Value = filter_th;
sldr_2.Value = dilate_rad;
% 
% hf = figure;
% set(hf, 'position', [500 10 1000 1000]);
% img_recall(t, im_ph_fst, xedge, yedge);
% txt_1 = uicontrol(hf,'Style', 'text',...
%         'String','Filter threshold',...
%         'Position', [250 120 100 20]);
% edt_1 = uicontrol(hf,'Style', 'edit',...
%         'String',num2str(filter_th),...
%         'Position', [250 80 100 20]);
% sldr_1 = uicontrol(hf,'Style','slider',...
%         'min',0.001,'max',1,'Value',filter_th, 'SliderStep',[0.0005 0.001],...
%         'callback','filter_th = sldr_1.Value; edt_1.String=filter_th; [~, xedge, yedge, domain] = find_boundary(im_ph_fst, filter_th, dilate_rad);img_recall(t, im_ph_fst, xedge, yedge);',...
%         'Position', [250 97 100 20]);
% txt_2 = uicontrol(hf,'Style', 'text',...
%         'String','Dilate radius',...
%         'Position', [650 120 100 20]);
% edt_2 = uicontrol(hf,'Style', 'edit',...
%         'String',num2str(dilate_rad),...
%          'Position', [650 80 100 20]);
% btn_dn = uicontrol(hf,'Style','pushbutton',...
%         'String','-',...
%         'callback','dilate_rad = dilate_rad-1;  edt_2.String=dilate_rad; [~, xedge, yedge, domain] = find_boundary(im_ph_fst, filter_th, dilate_rad); img_recall(t, im_ph_fst, xedge, yedge);',...
%         'position', [650 100 20 20]);
% btn_up = uicontrol(hf,'Style','pushbutton',...
%         'String','+',...
%         'callback','dilate_rad = dilate_rad+1; edt_2.String=dilate_rad; [~, xedge, yedge, domain] = find_boundary(im_ph_fst, filter_th, dilate_rad); img_recall(t, im_ph_fst, xedge, yedge);',...
%         'position', [730 100 20 20]);  
% btn_fn = uicontrol(hf,'Style','pushbutton',...
%         'String', 'Apply', ...
%         'Position', [800 80 50 20],...
%         'Callback', 'close;');
% waitfor(hf);

hf = figure;
set(hf, 'position', [500 10 1000 1000]);
img_recall(t, im_ph_fst, xedge, yedge);
txt_1 = uicontrol(hf,'Style', 'text',...
        'String','Filter threshold',...
        'Position', [250 120 100 20]);
edt_1 = uicontrol(hf,'Style', 'edit',...
        'String',num2str(filter_th),...
        'Position', [250 80 100 20]);
sldr_1 = uicontrol(hf,'Style','slider',...
        'min',0.001,'max',1,'Value',filter_th, 'SliderStep',[0.0005 0.001],...
        'callback','filter_th = sldr_1.Value; edt_1.String=filter_th; [~, xedge, yedge, domain] = find_boundary(im_ph_fst, filter_th, dilate_rad);img_recall(t, im_ph_fst, xedge, yedge);',...
        'Position', [250 97 100 20]);
txt_2 = uicontrol(hf,'Style', 'text',...
        'String','Dilate radius',...
        'Position', [650 120 100 20]);
edt_2 = uicontrol(hf,'Style', 'edit',...
        'String',num2str(dilate_rad),...
         'Position', [650 80 100 20]);
btn_dn = uicontrol(hf,'Style','pushbutton',...
        'String','-',...
        'callback','dilate_rad = dilate_rad-1;  edt_2.String=dilate_rad; [~, xedge, yedge, domain] = find_boundary(im_ph_fst, filter_th, dilate_rad); img_recall(t, im_ph_fst, xedge, yedge);',...
        'position', [650 100 20 20]);
btn_up = uicontrol(hf,'Style','pushbutton',...
        'String','+',...
        'callback','dilate_rad = dilate_rad+1; edt_2.String=dilate_rad; [~, xedge, yedge, domain] = find_boundary(im_ph_fst, filter_th, dilate_rad); img_recall(t, im_ph_fst, xedge, yedge);',...
        'position', [730 100 20 20]);  
btn_fn = uicontrol(hf,'Style','pushbutton',...
        'String', 'Apply', ...
        'Position', [800 80 50 20],...
        'Callback', 'close;');
btn_fn2 = uicontrol(hf,'Style','pushbutton',...
        'String', 'Modify', ...
        'Position', [800 120 50 20],...
        'Callback', 'close; modif = 1;');
btn_fn3 = uicontrol(hf,'Style','pushbutton',...
        'String', 'Draw', ...
        'Position', [800 160 50 20],...
        'Callback', 'close; modif = 2;');
waitfor(hf);
    if modif == 1
        hdmodif = figure;
        title(['1 Wait until line displaying / 2 Modify the line / 3 Close this window']);
        imshow(im_t); hold on
        j = boundary(xedge,yedge,0.5);
        h = impoly(gca, [xedge(j),yedge(j)]);
        setColor(h,'yellow');
        addNewPositionCallback(h,@(p) title(mat2str(p,3)));
        fcn = makeConstrainToRectFcn('impoly', get(gca,'XLim'), get(gca,'YLim'));
        setPositionConstraintFcn(h,fcn);
        position = wait(h);
        domain = poly2mask(position(:,1),position(:,2),M,N); 
        modif = 0;
        close;
    elseif modif == 2
        hdmodif = figure;
        title(['1 Wait until line displaying / 2 Draw the line / 3 Double click in the middle']);
        imshow(im_ph_fst); hold on
        position = wait(impoly());
        domain = poly2mask(position(:,1),position(:,2),M,N);
        modif = 0;
        close;
    end
    
[t_num, crpdimg_dir]=img_crop(t, t_num, img_dir, domain, im_ph_fst, M, N);
rmdir([img_dir '/Rawdata'], 's');
save('inp_setup.mat', 'th', 'E', 'nu', 'pix_size', 'grid_resolution_ph', 'grid_space_ph', 'overlap_ph','grid_resolution_bd', 'grid_space_bd', 'overlap_bd', 'int_time', 'img_dir', 'crpdimg_dir', 't_num', 'code_dir');

%% =========== 4_Calculating PIV of cell and bead displacement ============
PIV_calculation;

%% ====================== 5_Extrating domain images =====================
load([img_dir '/PIV_ph.mat']);

cd(crpdimg_dir);
clear btn_dn btn_fn btn_fn2 btn_fn3 btn_fn4 btn_up dhf edt_1 edt_2 sldr_1 txt_1 txt_2

while t >= 1 && t <= t_num
    disp(['======== ' num2str(t) '/' num2str(t_num) ' ========']);
    im_t = imread('phase_crpt.tif',t);
    [M, N] = size(im_t);
    [~, xedge, yedge, domain] = find_boundary(im_t, filter_th, dilate_rad);      
    dhf = figure;
    set(dhf, 'units','normalized','outerposition',[0 0 1 1]); %'position', [500 10 1000 1000]);
    img_recall(t, im_t, xedge, yedge);
    txt_1 = uicontrol(dhf,'Style', 'text',...
            'String','Filter threshold',...
            'Position', [250 120 100 20]);
    edt_1 = uicontrol(dhf,'Style', 'edit',...
            'String',num2str(filter_th),...
            'Position', [250 80 100 20]);
    sldr_1 = uicontrol(dhf,'Style','slider',...
            'min',0.001,'max',1,'Value',filter_th, 'SliderStep',[0.0005 0.001],...
            'callback','filter_th = sldr_1.Value; edt_1.String=filter_th; [~, xedge, yedge, domain] = find_boundary(im_t, filter_th, dilate_rad);img_recall(t, im_t, xedge, yedge);',...
            'Position', [250 97 100 20]);
    txt_2 = uicontrol(dhf,'Style', 'text',...
            'String','Dilate radius',...
            'Position', [650 120 100 20]);
    edt_2 = uicontrol(dhf,'Style', 'edit',...
            'String',num2str(dilate_rad),...
             'Position', [650 80 100 20]);
    btn_dn = uicontrol(dhf,'Style','pushbutton',...
            'String','-',...
            'callback','dilate_rad = dilate_rad-1;  edt_2.String=dilate_rad; [~, xedge, yedge, domain] = find_boundary(im_t, filter_th, dilate_rad); img_recall(t, im_t, xedge, yedge);',...
            'position', [650 100 20 20]);
    btn_up = uicontrol(dhf,'Style','pushbutton',...
            'String','+',...
            'callback','dilate_rad = dilate_rad+1; edt_2.String=dilate_rad; [~, xedge, yedge, domain] = find_boundary(im_t, filter_th, dilate_rad); img_recall(t, im_t, xedge, yedge);',...
            'position', [730 100 20 20]);   
    btn_fn = uicontrol(dhf,'Style','pushbutton',...
            'String', 'Apply', ...
            'Position', [800 80 50 20],...
            'Callback', 'close; modif = 0;');
    btn_fn2 = uicontrol(dhf,'Style','pushbutton',...
        'String', 'Modify', ...
        'Position', [800 120 50 20],...
        'Callback', 'close; modif = 1;');
    btn_fn3 = uicontrol(dhf,'Style','pushbutton',...
        'String', 'Draw', ...
        'Position', [800 160 50 20],...
        'Callback', 'close; modif = 2;');
    btn_fn4 = uicontrol(dhf,'Style','pushbutton',...
        'String', 'Back', ... 
        'Position', [800 200 50 20],...
        'Callback', 'close; modif = 0; t = t-1;');
    waitfor(dhf);
    
    if modif == 1
        hdmodif = figure('units','normalized','outerposition',[0 0 1 1]);
        ylabel(['Cell image (at ' num2str((t-1)*10) ' min)']);
        j = boundary(xedge,yedge,0.5);
        imshow(im_t,'Border','tight','InitialMagnification','fit'); hold on
        set(hdmodif,'menubar','none');
        addNewPositionCallback(h,@(p) title(mat2str(p,3)));
        fcn = makeConstrainToRectFcn('impoly', get(gca,'XLim'), get(gca,'YLim'));
        h = impoly(gca, [xedge(j),yedge(j)]);
        setColor(h,'yellow');
        setPositionConstraintFcn(h,fcn);
        position = wait(h);
        domain = poly2mask(position(:,1),position(:,2),M,N); 
        modif = 0;
        close;
    elseif modif == 2
        hdmodif = figure('units','normalized','outerposition',[0 0 1 1]);
        set(hdmodif,'menubar','none');hold on
        imshow(im_t,'Border','tight','InitialMagnification','fit'); hold on
        ylabel(['Cell image (at ' num2str((t-1)*10) ' min)']);
        position = wait(impoly());
        domain = poly2mask(position(:,1),position(:,2),M,N);
        modif = 0;
%         [domain,~,~] = impoly();
%         modif = 0;
        close;
    end
    
    domain_all(:,:,t) = domain;
    t = t+1;
    save('domain_backup.mat', 'domain_all', 't', 'pix_size', 'int_time', 'grid_resolution_ph', 'grid_space_ph', 'overlap_ph', 'img_dir', 'crpdimg_dir', 't_num', 'filter_th', 'dilate_rad');
    clear dhf txt_1 edt_1 sldr_1 txt_2 edt_2 btn_dn btn_up btn_fn btn_fn2 btn_fn3 btn_fn4 hdmodif;
end
%%  
% Disp_ph = sqrt(dx_ph.^2+dy_ph.^2);
% 
% if exist('t') ~= 1 || t == 1
%     t = 1;
%     domain_sm2 = false([size(x_ph,1) size(x_ph,2)]);
%     filter_th = 0.05; % Threshold for Matlab Sobel filter
%     dilate_rad = 7; % Radius of circular kernal used to dilate bright spots
% end
% 
% while t >= 1 && t <= t_num
%     disp(['======== ' num2str(t) '/' num2str(t_num) ' ========']);
%     im_t = imread('phase_crpt.tif',t);
%     [M, N] = size(im_t);
% 
%     if t ~= 1
%         domain_sm1 = imdilate(imresize(domain_all(:,:,t-1), [size(x_ph,1) size(x_ph,2)]),strel('disk',1,0));
%         [ph_grad, ~] = imgradient(dx_ph(:,:,t-1),dy_ph(:,:,t-1));
%         ph_edge = edge(ph_grad, 'Canny', 0.05);
%         ph_edge2 = imclose(ph_edge,strel('disk',2,0)) | domain_sm1;
%         ph_edge3 = imfill(ph_edge2,'holes');
%         domain_sm2 = ph_edge3;
%         CC = bwconncomp(ph_edge3);
%         numPixels = cellfun(@numel,CC.PixelIdxList);
%         [biggest,idx] = max(numPixels);
%         ph_edge3(CC.PixelIdxList{idx}) = 0;
%         domain_sm2 = domain_sm2 - ph_edge3;
% %         domain_sm = imdilate(imresize(domain, [size(x_ph,1) size(x_ph,2)]),strel('disk',1,0));
% %         B = bwboundaries(domain_all(:,:,1+t));
% %         [~, max_idx] = max(cellfun('size', B, 1));
% %         B_rowcol = flipud(B{max_idx});
% %         clear xedge_sm yedge_sm;
% %         xedge_sm = B_rowcol(:,2);
% %         yedge_sm = B_rowcol(:,1);
% %         xedge_sm = xedge_sm + griddata(x_ph,y_ph,dx_ph(:,:,t),xedge_sm,yedge_sm);
% %         yedge_sm = yedge_sm + griddata(x_ph,y_ph,dy_ph(:,:,t),xedge_sm,yedge_sm);
% %         domain_sm2 = inpolygon(x_ph,y_ph,xedge_sm,yedge_sm);
%     end
%     [xedge, yedge, domain] = find_boundary2(t, im_t, domain_sm2, filter_th, dilate_rad);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dhf = figure;
%     set(dhf, 'units','normalized','outerposition',[0 0 1 1]); %'position', [500 10 1000 1000]);
%     img_recall(t, im_t, xedge, yedge);
%     txt_1 = uicontrol(dhf,'Style', 'text',...
%             'String','Filter threshold',...
%             'Position', [250 120 100 20]);
%     edt_1 = uicontrol(dhf,'Style', 'edit',...
%             'String',num2str(filter_th),...
%             'Position', [250 80 100 20]);
%     sldr_1 = uicontrol(dhf,'Style','slider',...
%             'min',0.001,'max',1,'Value',filter_th, 'SliderStep',[0.0005 0.001],...
%             'callback','filter_th = sldr_1.Value; edt_1.String=filter_th; [xedge, yedge, domain] = find_boundary2(t, im_t, domain_sm2, filter_th, dilate_rad);img_recall(t, im_t, xedge, yedge);',...
%             'Position', [250 97 100 20]);
%     txt_2 = uicontrol(dhf,'Style', 'text',...
%             'String','Dilate radius',...
%             'Position', [650 120 100 20]);
%     edt_2 = uicontrol(dhf,'Style', 'edit',...
%             'String',num2str(dilate_rad),...
%              'Position', [650 80 100 20]);
%     btn_dn = uicontrol(dhf,'Style','pushbutton',...
%             'String','-',...
%             'callback','dilate_rad = dilate_rad-1;  edt_2.String=dilate_rad; [xedge, yedge, domain] = find_boundary2(t, im_t, domain_sm2, filter_th, dilate_rad); img_recall(t, im_t, xedge, yedge);',...
%             'position', [650 100 20 20]);
%     btn_up = uicontrol(dhf,'Style','pushbutton',...
%             'String','+',...
%             'callback','dilate_rad = dilate_rad+1; edt_2.String=dilate_rad; [xedge, yedge, domain] = find_boundary2(t, im_t, domain_sm2, filter_th, dilate_rad); img_recall(t, im_t, xedge, yedge);',...
%             'position', [730 100 20 20]);   
%     btn_fn = uicontrol(dhf,'Style','pushbutton',...
%             'String', 'Apply', ...
%             'Position', [800 80 50 20],...
%             'Callback', 'close; modif = 0;');
%     btn_fn2 = uicontrol(dhf,'Style','pushbutton',...
%         'String', 'Modify', ...
%         'Position', [800 120 50 20],...
%         'Callback', 'close; modif = 1;');
%     btn_fn3 = uicontrol(dhf,'Style','pushbutton',...
%         'String', 'Draw', ...
%         'Position', [800 160 50 20],...
%         'Callback', 'close; modif = 2;');
%     btn_fn4 = uicontrol(dhf,'Style','pushbutton',...
%         'String', 'Back', ... 
%         'Position', [800 200 50 20],...
%         'Callback', 'close; modif = 0; t = t-2;');
%     waitfor(dhf);
%     
%     if modif == 1
%         hdmodif = figure('units','normalized','outerposition',[0 0 1 1]);
%         imshow(im_t,'Border','tight','InitialMagnification','fit'); hold on
%         set(hdmodif,'menubar','none');
%         ylabel(['Cell image (at ' num2str((t-1)*10) ' min)']);
%         j = boundary(xedge,yedge,0.5);
%         h = impoly(gca, [xedge(j),yedge(j)]);
%         setColor(h,'yellow');
%         addNewPositionCallback(h,@(p) title(mat2str(p,3)));
%         fcn = makeConstrainToRectFcn('impoly', get(gca,'XLim'), get(gca,'YLim'));
%         setPositionConstraintFcn(h,fcn);
%         position = wait(h);
%         domain = poly2mask(position(:,1),position(:,2),M,N); 
%         modif = 0;
%         close;
%     elseif modif == 2
%         hdmodif = figure('units','normalized','outerposition',[0 0 1 1]);
%         set(hdmodif,'menubar','none');hold on
%         imshow(im_t,'Border','tight','InitialMagnification','fit'); hold on
%         ylabel(['Cell image (at ' num2str((t-1)*10) ' min)']);
%         position = wait(impoly());
%         domain = poly2mask(position(:,1),position(:,2),M,N);
%         modif = 0;
% %         [domain,~,~] = impoly();
% %         modif = 0;
%         close;
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     domain_all(:,:,t) = domain;
%     t = t+1;
%     save('domain_backup.mat', 'domain_all', 't', 'pix_size', 'int_time', 'grid_resolution_ph', 'grid_space_ph', 'overlap_ph', 'img_dir', 'crpdimg_dir', 't_num', 'filter_th', 'dilate_rad');
%     clear dhf txt_1 edt_1 sldr_1 txt_2 edt_2 btn_dn btn_up btn_fn btn_fn2 btn_fn3 btn_fn4 hdmodif;
% end
%%  
if exist('domain.tif','file') == 2
    delete('domain.tif'); 
end

if t < t_num
    t_num = t-1;
end

for t = 1:t_num
    imwrite(domain_all(:,:,t),'domain.tif','tif','writemode','append');
end

%% ============= 6_Extra inputs (editable point for 't_num') ==============
%%%%%%%%%%%%%%%%%%% for saving memory %%%%%%%%%%%%%%%%%%%

cd(img_dir);
save('inp_setup.mat', 'th', 'E', 'nu', 'pix_size', 'grid_resolution_ph', 'grid_space_ph', 'overlap_ph','grid_resolution_bd', 'grid_space_bd', 'overlap_bd', 'int_time', 'img_dir', 'crpdimg_dir', 't_num', 'code_dir');

clear;
close all;
clc;
load('inp_setup.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% =============== 7_Calculating traction force microscopy ================
traction;
 
%% ============== 8_Calculating monolayer stress microscopy ===============
%if exist('model.in','file') == 0
    copyfile([code_dir '\model.in'],[img_dir '/model.in']);
%end
 
stress;
 
%%% average normal stresses
load('traction.mat');
load('stresses.mat');
 
% stress_time in stresses.mat :  [sigma_xx    sigma_yy    sigma_xy    P_max    P_min    vx_max    vy_max    SE]
 
pixelTometer = pix_size*1.e-6; 
tn = zeros([size(x,1) size(x,2) t_num]);
sh = zeros([size(x,1) size(x,2) t_num]);
 
for t = 1:t_num
    st_max = stress_time{t}(:,4);
    st_min = stress_time{t}(:,5);

    st_x = round(node_coords_time{t}(:,1)/pixelTometer);
    st_y = round(node_coords_time{t}(:,2)/pixelTometer);

    avg_st = (st_max+st_min)/2;
    st_sh = (st_max-st_min)/2;
    %[angle_maxst, rho_maxst] = cart2pol();
    
    for stsize = 1:size(avg_st,1)
        [~, ii] = find(x==st_x(stsize),1);
        [jj, ~] = find(y==st_y(stsize),1);
        tn(jj, ii, t) = avg_st(stsize);
        sh(jj, ii, t) = st_sh(stsize);
    end
end
 
save('tension.mat','x', 'y', 'tn', 'sh');
 
%% ======================== 9_Radial coordination =========================
radial_coordination;
 
%% ========================== 10_Plotting maps =============================
plot_basic_modif;
 

%% ===================== 11_Calculating trajectories ======================
trajectory_edit;
 
%% ==================== 12_Autocorrelation of tension =====================
% autocorr_tension;
trajectory_edit2;
plot_basic_modif2
autocorr_tension2;
