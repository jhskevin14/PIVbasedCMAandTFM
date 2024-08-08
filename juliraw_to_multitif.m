function [img_dir]=juliraw_to_multitif
%% ============================= GUI setting ==============================

% assinging the folder of raw image files by Juli
disp('======== Select the folder of a point for analysis (Raw data by JuLI) ========');
Maindir = uigetdir('','Select the folder of a point for analysis (Raw data by JuLI)');
% assinging the folder of trypsin image files by Juli
disp('======== Select the bead trypsin image for analysis ========');
[Trimg, Trdir] = uigetfile('*.jpg','Select the bead trypsin image for analysis');
% bead color selection
BeadColor = questdlg('Which color of the fluorescent bead did you use?', 'Bead Color', 'GFP','RFP','GFP');
% assinging the folder for analysis
disp('======== Select or make a folder to save for analysis ========');
img_dir = uigetdir('','Select or make a folder to save for analysis');

%% ===================== Generating multi-tif images ======================

Phaseimg = dir(fullfile([Maindir '/Bright/*.jpg']));
Beadimg = dir(fullfile([Maindir '/' BeadColor '/*.jpg']));

t_num = numel(Phaseimg);

cd(img_dir);
mkdir('Rawdata');

Trimg = imread([Trdir Trimg]);
imwrite(Trimg(:,:,2),[img_dir '/Rawdata/bead_tr.tif'],'tif');

[M, N] = size(Trimg(:,:,2));
ph_t = uint8(zeros(M, N, t_num));
bd_t = uint8(zeros(M, N, t_num));

for t=1:t_num
        ph_0 = imread([Maindir '/BRIGHT/' Phaseimg(t).name]);
        ph_t(:,:,t) = ph_0(:,:,1);
        bd_0 = imread([Maindir '/' BeadColor '/' Beadimg(t).name]);
        bd_t(:,:,t) = bd_0(:,:,2);
        imwrite(ph_t(:,:,t),[img_dir '/Rawdata/phase.tif'],'tif','writemode','append');
        imwrite(bd_t(:,:,t),[img_dir '/Rawdata/bead.tif'],'tif','writemode','append');
end
end