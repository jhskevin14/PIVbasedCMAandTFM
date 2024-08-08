function autocorr_tension
clear;
close all;
clc;

load('inp_setup.mat');
load('tension.mat');
load('trajectories.mat');

if exist('half_tn','var') == 1
    clear atcorr_tn half_tn;
end

t_num = t_num-1;
node_num = size(traj_x,1);

for t = 1:t_num
    for node = 1:node_num
        [~,Y] = min(min(abs(x_cell-traj_x(node,t))));
        [~,X] = min(abs(y_cell-traj_y(node,t)));
        tension(node,t) = tn(X(1,1),Y(1,1),t);
    end
end

t_starting = 24; % starting point for sampling (at 240 min = after 4 hr)
t_sampling = 24; % sampling duration (400 min)
% atcorr_tn = zeros(node_num,t_num);

for node = 1:node_num
    atcorr_tn(node,:) = xcorr(tension(node,t_starting:t_starting+t_sampling), 'coeff');
%    atcorr_tn(node,:) = xcorr2(tension(node,t_starting:t_starting+t_sampling));
end

for node = 1:node_num
    if isnan(atcorr_tn(node,1))
        half_tn(node,1) = NaN;
    else
         half_tn(node,1) = find(atcorr_tn(node,:)>0.5,1,'last') - t_sampling;
%         [~,hp] = find(atcorr_tn(node,:) >= (max(atcorr_tn(node,:)))/2, 1, 'first');
%         half_tn(node,1) = t_sampling - hp;
    end
end
hftn = figure;
set(hftn, 'color','w','position', [500 200 640 640]);
scatter(traj_x(:,t_starting+t_sampling).*pix_size, traj_y(:,t_starting+t_sampling).*pix_size, 100, half_tn.*10, 'filled');
set(gcf,'PaperPositionMode','auto');
c = colorbar; c.Label.String = 'min'; c.Label.FontSize = 15;
set(gca,'CLim',[0 180]); 
colormap(hftn, jet);
title(['Correlation time of tension (' num2str(t_starting/6) ' hr ~ ' num2str((t_starting+t_sampling)/6) ' hr)'],'FontSize',15);
xlabel('\mum','fontsize',15);
ylabel('\mum','fontsize',15);
axis equal;
set(gca,'box','on','fontsize',11);
set(gca,'Ydir','reverse');
set(gca,'color','k');
axis ([-1000 1000 -1000 1000]);
% 
print([img_dir '/Figures/Autocorrelation_tension'], '-dpng','-r300','-r0');
close;
end