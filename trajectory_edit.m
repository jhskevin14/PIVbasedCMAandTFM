function trajectory_edit
% compute_cell_trajectories.m
%
% Compute pathlines for cell trajectories based on cell displacements
% calculated from image correlation.
%
% Algorithm adapted from a code written by Chan Young Park and Jacob Notbohm, Harvard School
% of Public Health, May 2012, and February 2014
% 
% Written by Jang, Hwanseok, Korea University, February 2017

clear;
close all;
clc;

code_dir = 'D:\HJang\Documents\Dropbox\MatlabCode_HS\161201_codeHS';
addpath(code_dir);

load('inp_setup.mat');
if exist('extra_inp.mat','file') == 2
    load('extra_inp.mat');
end
load('PIV_ph.mat');

t_num = t_num-1;

% Name of multipage tif file with domains
domain_path = [crpdimg_dir '/domain.tif'];

% Load data

% Rename variables associated with cell displacements so they aren't
% confused with substrate displacements
x_cell = x_ph;
y_cell = y_ph;
u_cell = dx_ph;
v_cell = dy_ph;

% Get center from center of first domain image
domain_fst = imread(domain_path,1);
domain_fst = double(domain_fst); % Convert to double precision
domain_fst = domain_fst/max(domain_fst(:)); % Set max value to 1

% Domain resize
[M, N] = size(x_ph);
domain_fst = imresize(domain_fst, [M N]);
domain_fst = round(domain_fst);

% Centroid coordinates
xc = sum(x_ph(:).*domain_fst(:)) / sum(domain_fst(:)); % Units: pix (same units as x and y)
yc = sum(y_ph(:).*domain_fst(:)) / sum(domain_fst(:));

% Convert to um
x_cell = (x_cell-xc)*pix_size;
y_cell = (y_cell-yc)*pix_size;
u_cell = u_cell*pix_size;
v_cell = v_cell*pix_size;


% Factor by which to downsample data. (Need to be careful here -- final
% number of grid points should match number of cells.)
fd = 2;

% Get initial grid points for trajectory computation
x_cell2=downsample(x_cell,fd); x_cell2=downsample(x_cell2',fd)';
y_cell2=downsample(y_cell,fd); y_cell2=downsample(y_cell2',fd)';

% First domain
domain_fst = imread(domain_path,1);

% Get domain's boundary coordaintes
BNDRY = bwboundaries(domain_fst); % This should be 1 cell with boundary coordinates
BNDRY = BNDRY{1};
x_bndry = BNDRY(:,2); % x coordinates are the columns
y_bndry = BNDRY(:,1); % y coordinates are the rows
% Shift at center and convert from pix to um
x_bndry = (x_bndry-xc)*pix_size; y_bndry = (y_bndry-yc)*pix_size;

% Get indices of points inside region given by (x_bndry,y_bndry)
IDX = inpolygon(x_cell2,y_cell2,x_bndry,y_bndry);

% Number of trajectories is given by number of nonzero elements in IDX
num_traj = sum(double(IDX(:)));

% Preallocate arrays for trajectories. These arrays start at each 
% downsampled  grid point within the domain. Rows are different traj and 
% columns are different times.
traj_x = nan*zeros(num_traj,t_num);
traj_y = nan*zeros(num_traj,t_num);

% First set of trajectories are given by points inside domain (ie, where IDX==1)
traj_x(:,1) = x_cell2(IDX);
traj_y(:,1) = y_cell2(IDX);

% Scan through remaining time points adding to the path coordinates
for t=2:t_num
    % k-th cell displacements
    u_cell_k = u_cell(:,:,t); % units: um
    v_cell_k = v_cell(:,:,t);
    
    % Interpolate k-th cell displacements to gridpoints from (k-1)th
    % timepoint
    displ_x = griddata(x_cell,y_cell,u_cell_k,traj_x(:,t-1),traj_y(:,t-1));
    displ_y = griddata(x_cell,y_cell,v_cell_k,traj_x(:,t-1),traj_y(:,t-1));
    
    % Add to trajectory arrays
    traj_x(:,t) = traj_x(:,t-1) + displ_x; % Units: um
    traj_y(:,t) = traj_y(:,t-1) + displ_y;
end

% Save data
save('trajectories.mat','x_cell','y_cell','traj_x','traj_y');
%% Compute Distance / Path length ratio

% Number of trajectories, N, and number of time points, K
[N, K] = size(traj_x);

for t=2:K
    pathlength(:,t-1) = sqrt((traj_x(:,t)-traj_x(:,t-1)).^2+(traj_y(:,t)-traj_y(:,t-1)).^2);
end
for p=1:size(pathlength,1)
    pathlength_2(p) = sum(pathlength(p,:));
end
distance = sqrt((traj_x(:,K)-traj_x(:,1)).^2+(traj_y(:,K)-traj_y(:,1)).^2);
distance = transpose(distance);
dp_ratio = distance./pathlength_2;

%% Compute D2min

% Distance between cells n and j to use in calculation of D2min
% D2min_rad = 40; % Units: um
% Instead of distance between cells n and j, specify number of nearest
% neighbors. (Choosing number of nearest neighbors is better, because the
% cells get farther away from each other over time.)
n_nearestneighbors = 8; % Nader suggests 8
% Limit for color plots
D2min_lim = 1000; % Units: um^2
% Time shift used to compute Delta d_{ij} in equation above. Usually 
% pick the time required to move a distance of 1 cell diameter.
deltaT = 10; % (for 100 min) Units: time increments
% Time between images
time_increment = 10; % Units: min



% Preallocate matrix of values
D2min = zeros(N,K-deltaT)*nan;
for k=1:K-deltaT
    % Trajectory positions at time point k
    traj_x_k = traj_x(:,k);
    traj_y_k = traj_y(:,k);
    % Trajectory positions at time point k+1. Note that
    % D2min is often computed over a period of time equal to the time
    % required for a cell to move one cell diameter. Call this time deltaT,
    % where deltaT is a specific number of time steps
    traj_x_kplus = traj_x(:,k+deltaT);
    traj_y_kplus = traj_y(:,k+deltaT);
    
    % Compute D2min for each trajectory
    for i=1:N % i indexes over cell trajectories
        % Coordinates of cell i
        xi = traj_x_k(i); % Units: um
        yi = traj_y_k(i);
        
        xi_plus = traj_x_kplus(i);
        yi_plus = traj_y_kplus(i);
        
        % Get cells within a radius D2min_rad of (xi,yi)
        dist = sqrt((traj_x_k-xi).^2 + (traj_y_k-yi).^2); % Dist is a vector
        %         idx = dist<D2min_rad & dist>0;
        
        % Instead of using a distance, specify number of nearest neighbors
        [~, idx] = sort(dist,'ascend'); % Sort distances
        idx = idx(1:n_nearestneighbors); % Get n nearest neighbors
        
        xj = traj_x_k(idx); % Units: um
        yj = traj_y_k(idx);
        
        xj_plus = traj_x_kplus(idx);
        yj_plus = traj_y_kplus(idx);
        
        % Compute Delta d_ij
        Delta_dx = (xj_plus-xi_plus) - (xj-xi);
        Delta_dy = (yj_plus-yi_plus) - (yj-yi);
        
        % Compute strain tensor E that best matches the motion of cell i
        % and the cells indexed by j. I compute E from Eqs. 2.12-2.14 of
        % Falk & Langer, PRE 57, 6, 1998. 
        % Double check that these equations make sense!!
        X(1,1) = sum((xj_plus-xi_plus).*(xj-xi));
        X(1,2) = sum((xj_plus-xi_plus).*(yj-yi));
        X(2,1) = sum((yj_plus-yi_plus).*(xj-xi));
        X(2,2) = sum((yj_plus-yi_plus).*(yj-yi));
        
        Y(1,1) = sum((xj-xi).*(xj-xi));
        Y(1,2) = sum((xj-xi).*(yj-yi));
        Y(2,1) = sum((yj-yi).*(xj-xi));
        Y(2,2) = sum((yj-yi).*(yj-yi));
        
        E = X*inv(Y) - eye(2); % * should be matrix multiplication; eye(2) is 2x2 identity matrix
        
        % Compute "affine displacement" by computing strain E times d_ij
        Edx = E(1,1)*(xj-xi) + E(1,2)*(yj-yi);
        Edy = E(2,1)*(xj-xi) + E(2,2)*(yj-yi);
        
        % Take the mean over all indices j to get D2min(i)
        D2min(i,k) = mean( (Delta_dx-Edx).^2 + (Delta_dy-Edy).^2 );       
    end
end

D2min_mean = mean(D2min,2);

%% --- PLOT TRAJECTORIES and Distance/Path length Ratio and D2min---

hf = figure;
set(hf,'color','w','position', [0 0 1920 850]);
set(gcf,'PaperPositionMode','auto');

fg1 = subplot(1,3,1);
pathlength_3 = zeros([size(pathlength,1) 1]);
for K=1:t-1
    pathlength_3 = pathlength_3 + pathlength(:,K);
    scatter(traj_x(:,K), traj_y(:,K), 2, pathlength_3, 'filled');
    hold on;
end
hold off;
c = colorbar; c.Label.String = 'Path length (\mum)'; c.Label.FontSize = 15;
set(gca,'CLim',[0 300]); colormap(fg1, cool);
title('Trajectories with path length','FontSize',15);
xlabel('\mum','fontsize',15);
ylabel('\mum','fontsize',15);
axis equal;
set(gca,'box','on','fontsize',11);
set(gca,'Ydir','reverse');
axis ([-800 800 -800 800]);

fg2 = subplot(1,3,2);
scatter(traj_x(:,K), traj_y(:,K), 50, dp_ratio, 'filled')
c = colorbar; c.Label.String = 'Distance/Path length ratio'; c.Label.FontSize = 15;
set(gca,'CLim',[0 1]); colormap(fg2, hot);
title('Distance/Path length','FontSize',15);
xlabel('\mum','fontsize',15);
ylabel('\mum','fontsize',15);
axis equal
set(gca,'box','on','fontsize',11);
set(gca,'Ydir','reverse');
axis ([-800 800 -800 800]);

fg3 = subplot(1,3,3);
scatter(traj_x(:,K), traj_y(:,K), 50, D2min_mean, 'filled')
c = colorbar; c.Label.String = '\mum^2'; c.Label.FontSize = 15;
set(gca,'CLim',[0 400]); 
colormap(fg3, jet);
title('Deviation from average local velocity (D^2_{min})','FontSize',15);
xlabel('\mum','fontsize',15);
ylabel('\mum','fontsize',15);
axis equal
set(gca,'box','on','fontsize',11);
set(gca,'Ydir','reverse');
axis ([-800 800 -800 800]);

print('-dpng', '-r0', [img_dir '/Figures/Traje_Dist-Pathlength']);
close;

end





