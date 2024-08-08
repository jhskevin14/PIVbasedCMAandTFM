function radial_coordination

clear;
close all;
clc;

load('inp_setup.mat');
if exist('extra_inp.mat','file') == 2
    load('extra_inp.mat');
end
load('PIV_ph.mat');
load('traction.mat');
load('tension.mat');

dx_ph = dx_ph - median(dx_ph(~isnan(dx_ph)));
dy_ph = dy_ph - median(dy_ph(~isnan(dy_ph)));

[M_vel, N_vel, Vel_rad, Vel_rad_aver, Vel_Theta, Vel_rho] = radial(x_ph, y_ph, dx_ph, dy_ph);
[M_tr, N_tr, Tr_rad, Tr_rad_aver, Tr_Theta, Tr_rho] = radial(x, y, tx, ty);
[M_tn, N_tn, Tn_rad, Tn_rad_aver, Tn_Theta, Tn_rho] = radial(x, y, tn, tn);

Vel_rad = Vel_rad*pix_size/10;
Vel_rad_aver = Vel_rad_aver*pix_size/10;
Vel_rho = Vel_rho*pix_size/10;

Tn_rad = Tn_rad/sqrt(2);
Tn_rad_aver = Tn_rad_aver/sqrt(2);
Tn_rho = Tn_rho/sqrt(2);

save('radial_coordination.mat','M_vel','N_vel', 'Vel_rad', 'Vel_rad_aver', 'Vel_Theta', 'Vel_rho',...
                               'M_tr','N_tr','Tr_rad', 'Tr_rad_aver', 'Tr_Theta', 'Tr_rho',...
                               'M_tn','N_tn','Tn_rad', 'Tn_rad_aver', 'Tn_Theta', 'Tn_rho');
end
function [M, N, Vec_rad, Vec_rad_aver, Theta, Rho] = radial(X, Y, DX, DY)

load('inp_setup.mat');
if exist('extra_inp.mat','file') == 2
    load('extra_inp.mat');
end

for t = 1:t_num
    domain_all(:,:,t) = imread([crpdimg_dir '/domain.tif'],t);    
end

[M, N, t_n] = size(DX);
if t_n > t_num
    t_n = t_num;
end
if exist('domain_resize','var') == 1
    clear domain_resize Theta Rho;
end
domain_resize(:,:,:) = imresize(domain_all(:,:,:), [M, N]);

s = regionprops(domain_resize(:,:,1),'centroid', 'Area');

[~,index] = max([s.Area]);
cent = cat(1, s(index).Centroid);
cent_r = round(cent);

xc = cent_r(1);
yc = cent_r(2);

Theta_glov = atan2(Y - Y(yc,xc), X - X(yc,xc));
Theta_glov = repmat(Theta_glov, [1 1 t_n]);

DX = DX(:,:,1:t_n).*domain_resize(:,:,1:t_n);
DY = DY(:,:,1:t_n).*domain_resize(:,:,1:t_n);

[Theta, Rho] = cart2pol(DX,DY);
Vec_rad = Rho.*cos(Theta - Theta_glov);

distanceToUL = sqrt((1-yc)^2 + (1-xc)^2);
distanceToUR = sqrt((1-yc)^2 + (M-xc)^2);
distanceToLL = sqrt((N-yc)^2 + (1-xc)^2);
distanceToLR= sqrt((N-yc)^2 + (M-xc)^2);
maxDistance = ceil(max([distanceToUL, distanceToUR, distanceToLL, distanceToLR]));

Vec_rad_Sums = zeros(maxDistance, t_n);
Vec_rad_Counts = zeros(maxDistance, t_n);

for t = 1 : t_n
    for i = 1 : N
        for j = 1 : M
            thisDistance = round(sqrt((j-yc)^2 + (i-xc)^2));
            if thisDistance <= 0
                continue;
            end
            Vec_rad_Sums(thisDistance,t) = Vec_rad_Sums(thisDistance, t) + Vec_rad(j,i,t);
            Vec_rad_Counts(thisDistance, t) = Vec_rad_Counts(thisDistance, t) + 1;
        end
    end
end

Vec_rad_aver = Vec_rad_Sums ./ Vec_rad_Counts;
end