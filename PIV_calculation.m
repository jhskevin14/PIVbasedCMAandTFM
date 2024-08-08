function PIV_calculation
%
load('inp_setup.mat');

%% ========================== Calling the images ==========================

cd(crpdimg_dir);
im_tr = imread('bead_tr_crpt.tif');
[M, N] = size(im_tr);
im_ph = uint8(zeros(M, N, t_num));
im_bd = uint8(zeros(M, N, t_num));

for t = 1:t_num
    im_ph(:,:,t) = imread('phase_crpt.tif',t);
    im_bd(:,:,t) = imread('bead_crpt.tif',t);
end

%% ============ Calculating migration velocity by phase images ============
disp('======== Starting calculation of deformation for migration velocity ========');

x_ph = cell(1,t_num-1);
y_ph = x_ph;
dx_ph = x_ph;
dy_ph = x_ph;

for t = 1:t_num-1
    tic;
    [x, y, dx, dy, ~] = disp_on_blocks_fast(double(im_ph(:,:,t)), ...   % reference image
                          double(im_ph(:,:,t+1)), ...                         % current image
                          grid_resolution_ph, ...                           % image size
                          grid_resolution_ph, ...                           % image size
                          overlap_ph, ...                                   % overlap between blocks
                          0, ...                                            % padding for the PIV, usually 0
                          'hanning', ...                                    % window
                          '2DCentroid', ...                                 % Subpixel interpolation method
                          0, ...                                            % threshold of center of mass calculation
                          2, ...                                           % Size of the window
                          4);                                               % number of iterations
    x_ph{t} = x;
    y_ph{t} = y;
    dx_ph{t} = dx;
    dy_ph{t} = dy;
    
    tm = toc;
    disp(['[Cell displacement] Correlation number ',num2str(t),' completed in ',num2str(tm/60),' mins.'])
end

x_ph_a(:,:) = x_ph{1};
x_ph = x_ph_a;
y_ph_a(:,:) = y_ph{1};
y_ph = y_ph_a;

for i=1:t
    dx_ph_a(:,:,i) = dx_ph{i};
    dy_ph_a(:,:,i) = dy_ph{i};
end

dx_ph = dx_ph_a;
dy_ph = dy_ph_a;

cd(img_dir);
save('PIV_ph.mat','x_ph','y_ph','dx_ph','dy_ph');

%% ============= Calculating bead displacement by bead images =============
disp('======== Starting calculation of deformation for traction calculation ========');

x_bd = cell(1,t_num);
y_bd = x_bd;
dx_bd = x_bd;
dy_bd = x_bd;


for t = 1:t_num
    tic;
    [x, y, dx, dy, ~] = disp_on_blocks_fast(double(im_tr(:,:)), ...   % reference image
                          double(im_bd(:,:,t)), ...                         % current image
                          grid_resolution_bd, ...                           % image size
                          grid_resolution_bd, ...                           % image size
                          overlap_bd, ...                                   % overlap between blocks
                          0, ...                                            % padding for the PIV, usually 0
                          'hanning', ...                                    % window
                          '2DCentroid', ...                                 % Subpixel interpolation method
                          0, ...                                            % threshold of center of mass calculation
                          2, ...                                           % Size of the window
                          4); 
    x_bd{t} = x;
    y_bd{t} = y;
    dx_bd{t} = dx;
    dy_bd{t} = dy;
    
    tm = toc;
    disp(['[Bead displacement] Correlation number ',num2str(t),' completed in ',num2str(tm/60),' mins.'])
end

x_bd_a(:,:) = x_bd{1};
x_bd = x_bd_a;
y_bd_a(:,:) = y_bd{1};
y_bd = y_bd_a;

for i=1:t
    dx_bd_a(:,:,i) = dx_bd{i};
    dy_bd_a(:,:,i) = dy_bd{i};
end

dx_bd = dx_bd_a;
dy_bd = dy_bd_a;

save('PIV_bd.mat','x_bd','y_bd','dx_bd','dy_bd');

end