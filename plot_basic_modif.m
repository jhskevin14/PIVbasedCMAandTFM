function plot_basic_modif

clear;
close all;
clc;

if exist('Figures', 'dir') == 0
    mkdir('Figures');
end
load('inp_setup.mat');
if exist('extra_inp.mat','file') == 2
    load('extra_inp.mat');
end
load('PIV_ph.mat');
load('traction.mat');
load('radial_coordination');
load('stresses.mat');
load('tension.mat');

%%% domain
for t = 1:t_num
    domain_all(:,:,t) = imread([crpdimg_dir '/domain.tif'],t);    
end
%%%

[M_ph, N_ph, t_ph] = size(dx_ph);
[M_bd, N_bd, t_bd] = size(dx);
domain_resize_ph(:,:,:) = imresize(domain_all(:,:,:), [M_ph, N_ph]);
domain_resize_bd(:,:,:) = imresize(domain_all(:,:,:), [M_bd, N_bd]);

t_max = t_num-1;
%% =================== RED-WHITE-BLUE colormap ===================
map = [repmat([0,0,1],[32,1]);[0,0,0];repmat([1,0,0],[32,1])];
r_ph = repmat(abs(linspace(1,-1,65)).^1,[3,1])';
map = map.*r_ph + 1 - r_ph;
%% =================== RED-Black-BLUE colormap ===================
load('map2.mat');
%% ======================= Displacement =========================
dx = dx - median(dx(~isnan(dx)));
dy = dy - median(dy(~isnan(dy)));
dx_ph = dx_ph - median(dx_ph(~isnan(dx_ph)));
dy_ph = dy_ph - median(dy_ph(~isnan(dy_ph)));

dx = dx(:,:,1:t_max).*domain_resize_bd(:,:,1:t_max);
dy = dy(:,:,1:t_max).*domain_resize_bd(:,:,1:t_max);
dx_ph = dx_ph(:,:,1:t_max).*domain_resize_ph(:,:,1:t_max);
dy_ph = dy_ph(:,:,1:t_max).*domain_resize_ph(:,:,1:t_max);
    
xq = imresize(x_ph,[20 20]);
yq = imresize(y_ph,[20 20]);

DispXq_bd = imresize(dx,[20 20]);
DispYq_bd = imresize(dy,[20 20]);
Disp_bd = sqrt(dx.^2+dy.^2);
dispLim_bd = max(max(Disp_bd));

DispXq_ph = imresize(dx_ph,[20 20]);
DispYq_ph = imresize(dy_ph,[20 20]);
Disp_ph = sqrt(dx_ph.^2+dy_ph.^2);
dispLim_ph = max(max(Disp_ph))*0.5;

x1_bd = dx(:,:,1:t_max);
y1_bd = dy(:,:,1:t_max);
x2_ph = dx_ph(:,:,1:t_max)*pix_size;
y2_ph = dy_ph(:,:,1:t_max)*pix_size;
theta_ellips = acosd((x1_bd.*x2_ph + y1_bd.*y2_ph)./(sqrt(x1_bd.^2+y1_bd.^2).*sqrt(x2_ph.^2+y2_ph.^2)));
theta_rad = acos((x1_bd.*x2_ph + y1_bd.*y2_ph)./(sqrt(x1_bd.^2+y1_bd.^2).*sqrt(x2_ph.^2+y2_ph.^2)));

sp = fix(t_max/5);

%% Cell_displacements
a=0; b=0; c=0;
fig_CellDisp = figure;
set(fig_CellDisp,'color','w','position', [55 20 1500 1125])
set(gcf,'PaperPositionMode','auto');
for t=1:sp:t_max
    a=a+1;
    if a == 4
        b=0; c=1;
    end
    sub_plot = subplot(2,3,a);
    set(sub_plot, 'position', [0.025+(0.325*b) 0.55-(0.5*c) 0.3 0.4]);
    pcolor(x_ph,y_ph,Disp_ph(:,:,t)*pix_size);
    colormap jet; shading flat;
    axis image; axis off; 
    set(gca,'Ydir','reverse');
    set(gca,'CLim',[0 10]); colorbar;
    hold on;
    hq = quiver(xq, yq, DispXq_ph(:,:,t)*pix_size,DispYq_ph(:,:,t)*pix_size);
    set(hq,'linewidth',1, 'Color', 'k');
    if t_num < 30
        title({'Cell displacement', ['(\mum, ' num2str((t-1)*10) 'min)']});
    elseif t_num >= 30
        title({'Cell displacement', ['(\mum, ' num2str(round((t-1)/6)) ' hours)']});
    end
    b = b+1;
    if a == 6
        break
    end
end
print([img_dir '/Figures/Cell_displacements'], '-dpng','-r300','-r0');
close;
%% Average speed in radial coordination
a=0; b=0; c=0;
fig_CellSpeedR = figure;
set(fig_CellSpeedR,'color','w','position', [55 20 1500 1125])
set(gcf,'PaperPositionMode','auto');
for t=1:sp:t_max
    a=a+1;
    if a == 4
        b=0; c=1;
    end
    sub_plot = subplot(2,3,a);
    set(sub_plot, 'position', [0.025+(0.325*b) 0.55-(0.5*c) 0.3 0.4]);
    pcolor(x_ph,y_ph,Vel_rad(:,:,t));
    colormap(map); shading flat;
    axis image; axis off; colorbar; 
    set(gca,'Ydir','reverse');
    set(gca,'CLim',[-0.8 0.8]); %set(gca,'CLim',[-0.5 0.5]);
    hold on;
    hq = quiver(xq, yq, DispXq_ph(:,:,t)*pix_size,DispYq_ph(:,:,t)*pix_size);
    set(hq,'linewidth',1, 'Color', 'k');
    if t_num < 30
        title({'Cellular velocity_{radial}', ['(\mum/min, ' num2str((t-1)*10) 'min)']});
    elseif t_num >= 30
        title({'Cellular velocity_{radial}', ['(\mum/min, ' num2str(round((t-1)/6)) ' hours)']});
    end
    b = b+1;
    if a == 6
        break
    end
end
print([img_dir '/Figures/Cell_velocity_radial'], '-dpng','-r300','-r0');
close;
%% Bead_displacements
a=0; b=0; c=0;
fig_BeadDisp = figure;
set(fig_BeadDisp,'color','w','position', [55 20 1500 1125]);
set(gcf,'PaperPositionMode','auto');
for t=1:sp:t_max
    a=a+1;
    if a == 4
        b=0; c=1;
    end
    sub_plot = subplot(2,3,a);
    set(sub_plot, 'position', [0.025+(0.325*b) 0.55-(0.5*c) 0.3 0.4]);
    pcolor(x,y,Disp_bd(:,:,t));
    colormap jet; shading flat;
    axis image; axis off; 
    set(gca,'Ydir','reverse');
    set(gca,'CLim',[0 0.5]); colorbar; %dispLim_bd(1,1,1)
    hold on;
    hq = quiver(xq, yq, DispXq_bd(:,:,t),DispYq_bd(:,:,t));
    set(hq,'linewidth',1, 'Color', 'k');
    if t_num < 30
        title({'Bead displacement', ['(\mu, ' num2str((t-1)*10) 'min)']});
    elseif t_num >= 30
        title({'Bead displacement', ['(\mu, '  num2str(round((t-1)/6)) ' hours)']});
    end
    b = b+1;
    if a == 6
        break
    end
end
print([img_dir '/Figures/Bead_displacements'], '-dpng','-r300','-r0');
close;

%% Traction_(magitude)
a=0; b=0; c=0;

% maxtrall = ceil(mean(Tr_rho(~(Tr_rho == 0)))*0.1)*20;

fig_Traction = figure;
set(fig_Traction,'color','w','position', [55 20 1500 1125]);
set(gcf,'PaperPositionMode','auto');
for t=1:sp:t_max
    a=a+1;
    if a == 4
        b=0; c=1;
    end
    sub_plot = subplot(2,3,a);
    set(sub_plot, 'position', [0.025+(0.325*b) 0.55-(0.5*c) 0.3 0.4]);
    pcolor(x, y, Tr_rho(:,:,t)); 
    colormap hot; shading interp;
    axis image; axis off; 
    set(gca,'Ydir','reverse');
    set(gca,'CLim',[0 50]); colorbar;
    hold on;
    hv = quiver(xq, yq, DispXq_bd(:,:,t),DispYq_bd(:,:,t));
    set(hv,'linewidth',1, 'Color', 'w');
    if t_num < 30
        title({'Traction force magnitude', ['(Pa, ' num2str((t-1)*10) 'min)']});
    elseif t_num >= 30
        title({'Traction force magnitude', ['(Pa, ' num2str(round((t-1)/6)) ' hours)']});
    end
    b = b+1;
    if a == 6
        break
    end
end
print([img_dir '/Figures/Traction_magnitude'], '-dpng','-r300','-r0');
close;

%% Traction_(Y-axis)
a=0; b=0; c=0;

% maxtrall = ceil(mean(Tr_rho(~(Tr_rho == 0)))*0.1)*20;

fig_Traction = figure;
set(fig_Traction,'color','w','position', [55 20 1500 1125]);
set(gcf,'PaperPositionMode','auto');
for t=1:sp:t_max
    a=a+1;
    if a == 4
        b=0; c=1;
    end
    sub_plot = subplot(2,3,a);
    set(sub_plot, 'position', [0.025+(0.325*b) 0.55-(0.5*c) 0.3 0.4]);
    pcolor(x, y, ty(:,:,t)*(-1)); 
    colormap(map2); shading interp;
    axis image; axis off; 
    set(gca,'Ydir','reverse');
    set(gca,'CLim',[-50 50]); colorbar;
    hold on;
    hv = quiver(xq, yq, DispXq_bd(:,:,t),DispYq_bd(:,:,t));
    set(hv,'linewidth',1, 'Color', 'k');
    if t_num < 30
        title({'Traction force_{Y-axis}', ['(Pa, ' num2str((t-1)*10) 'min)']});
    elseif t_num >= 30
        title({'Traction force_{Y-axis}', ['(Pa, ' num2str(round((t-1)/6)) ' hours)']});
    end
    b = b+1;
    if a == 6
        break
    end
end
print([img_dir '/Figures/Traction_Yaxis'], '-dpng','-r300','-r0');
close;

%% Traction_(Radial)
a=0; b=0; c=0;
% maxtrall = round(mean(Tr_rho(~(Tr_rad==0)))*0.1)*45;
fig_TractionR = figure;
set(fig_TractionR,'color','w','position', [55 20 1500 1125]);
set(gcf,'PaperPositionMode','auto');
for t=1:sp:t_max
    a=a+1;
    if a == 4
        b=0; c=1;
    end
    sub_plot = subplot(2,3,a);
    set(sub_plot, 'position', [0.025+(0.325*b) 0.55-(0.5*c) 0.3 0.4]);
    pcolor(x, y, Tr_rad(:,:,t)); 
    colormap(map2); shading interp;
    axis image; axis off; colorbar; 
    set(gca,'Ydir','reverse');
    set(gca,'CLim',[-50 50]);  % set(gca,'CLim',[-150 150]);
    hold on;
    hv = quiver(xq, yq, DispXq_bd(:,:,t),DispYq_bd(:,:,t));
    set(hv,'linewidth',1, 'Color', 'w');
    if t_num < 30
        title({'Traction force magnitude', ['(Pa, ' num2str((t-1)*10) 'min)']});
    elseif t_num >= 30
        title({'Traction force magnitude', ['(Pa, ' num2str(round((t-1)/6)) ' hours)']});
    end
    b = b+1;
    if a == 6
        break
    end
end
print([img_dir '/Figures/Traction_radial'], '-dpng','-r300','-r0');
close;
%% Tension (Average normal stress)
a=0; b=0; c=0;
fig_Tension = figure;
set(fig_Tension,'color','w','position', [55 20 1500 1125]);
set(gcf,'PaperPositionMode','auto');

for t=1:sp:t_max
    a=a+1;
    if a == 4
        b=0; c=1;
    end
    sub_plot = subplot(2,3,a);
    set(sub_plot, 'position', [0.025+(0.325*b) 0.55-(0.5*c) 0.3 0.4]);
%     surf(x, y, Tn_rho(:,:,t), 'edgecolor','none' );
%     view(2); colormap hot; shading interp;
    pcolor(x, y, Tn_rho(:,:,t));
    colormap hot; shading interp;
    % xlabel('\mum'); ylabel('\mum');
    axis image; axis off;
    %c = colorbar; c.Label.String = 'Pa';
    colorbar;
    set(gca,'CLim',[0 500]);
    set(gca,'Ydir','reverse');
    hold on;
    hv = quiver(xq, yq, DispXq_bd(:,:,t),DispYq_bd(:,:,t));
    set(hv,'linewidth',1, 'Color', 'k');
    hq = quiver(xq, yq, DispXq_ph(:,:,t)*pix_size,DispYq_ph(:,:,t)*pix_size);
    set(hq,'linewidth',1, 'Color', 'b');
    if t_num < 30
        title({'Tensional stress', ['(Pa, ' num2str((t-1)*10) 'min)']});
    elseif t_num >= 30
        title({'Tensional stress', ['(Pa, ' num2str(round((t-1)/6)) ' hours)']});
    end
    b = b+1;
    if a == 6
        break
    end
end
print([img_dir '/Figures/Tension'], '-dpng','-r300','-r0');
close;

%% Max shear stress
a=0; b=0; c=0;
fig_Shear = figure;
set(fig_Shear,'color','w','position', [55 20 1500 1125]);
set(gcf,'PaperPositionMode','auto');

for t=1:sp:t_max
    a=a+1;
    if a == 4
        b=0; c=1;
    end
    sub_plot = subplot(2,3,a);
    set(sub_plot, 'position', [0.025+(0.325*b) 0.55-(0.5*c) 0.3 0.4]);
%     surf(x, y, Tn_rho(:,:,t), 'edgecolor','none' );
%     view(2); colormap hot; shading interp;
    pcolor(x, y, sh(:,:,t));
    colormap hot; shading interp;
    % xlabel('\mum'); ylabel('\mum');
    axis image; axis off;
    %c = colorbar; c.Label.String = 'Pa';
    colorbar;
    set(gca,'CLim',[0 200]);
    set(gca,'Ydir','reverse');
    hold on;
    hv = quiver(xq, yq, DispXq_bd(:,:,t),DispYq_bd(:,:,t));
    set(hv,'linewidth',1, 'Color', 'k');
    hq = quiver(xq, yq, DispXq_ph(:,:,t)*pix_size,DispYq_ph(:,:,t)*pix_size);
    set(hq,'linewidth',1, 'Color', 'b');
    if t_num < 30
        title({'Max shear stress', ['(Pa, ' num2str((t-1)*10) 'min)']});
    elseif t_num >= 30
        title({'Max shear stress', ['(Pa, ' num2str(round((t-1)/6)) ' hours)']});
    end
    b = b+1;
    if a == 6
        break
    end
end
print([img_dir '/Figures/Max_shear'], '-dpng','-r300','-r0');
close;

%% Tension / Max shear stress
a=0; b=0; c=0;
fig_TnSh = figure;
set(fig_TnSh,'color','w','position', [55 20 1500 1125]);
set(gcf,'PaperPositionMode','auto');

for t=1:sp:t_max
    a=a+1;
    if a == 4
        b=0; c=1;
    end
    sub_plot = subplot(2,3,a);
    set(sub_plot, 'position', [0.025+(0.325*b) 0.55-(0.5*c) 0.3 0.4]);
%     surf(x, y, Tn_rho(:,:,t), 'edgecolor','none' );
%     view(2); colormap hot; shading interp;
    pcolor(x, y, Tn_rho(:,:,t)./sh(:,:,t));
    colormap jet; shading interp;
    % xlabel('\mum'); ylabel('\mum');
    axis image; axis off;
    %c = colorbar; c.Label.String = 'Pa';
    colorbar;
    set(gca,'CLim',[0 20]);
    set(gca,'Ydir','reverse');
    hold on;
    hv = quiver(xq, yq, DispXq_bd(:,:,t),DispYq_bd(:,:,t));
    set(hv,'linewidth',1, 'Color', 'k');
    hq = quiver(xq, yq, DispXq_ph(:,:,t)*pix_size,DispYq_ph(:,:,t)*pix_size);
    set(hq,'linewidth',1, 'Color', 'b');
    if t_num < 30
        title({'Tension/Max shear stress', ['(R, ' num2str((t-1)*10) 'min)']});
    elseif t_num >= 30
        title({'Tension/Max shear stress', ['(R, ' num2str(round((t-1)/6)) ' hours)']});
    end
    b = b+1;
    if a == 6
        break
    end
end
print([img_dir '/Figures/Tn_Sh'], '-dpng','-r300','-r0');
close;

%% Vel_trac_angle
a=0; b=0; c=0;
fig_vel_trac_angle = figure;
set(fig_vel_trac_angle,'color','w','position', [55 20 1500 1125]);
set(gcf,'PaperPositionMode','auto');
% 
% x1_bd = DispXq_bd(:,:,1:t_max);
% y1_bd = DispYq_bd(:,:,1:t_max);
% x2_ph = DispXq_ph(:,:,1:t_max)*pix_size;
% y2_ph = DispYq_ph(:,:,1:t_max)*pix_size;

% 
% a_ellip = 1/2*sqrt((x2_ph-x1_bd).^2+(y2_ph-y1_bd).^2);
% b_ellip = a_ellip.*tan(theta_ellips);
% 
% % lsp = linspace(0,2*pi);
% X_ellip = a_ellip.*cos(theta_ellips);
% Y_ellip = b_ellip.*sin(lsp);
% w_ellip = atan2(y2_ph-y1_bd,x2_ph-x1_bd);
% x_ellip = (x1_bd+x2_ph)./2 + X_ellip.*cos(w_ellip) - Y_ellip.*sin(w_ellip);
% y_ellip = (y1_bd+y2_ph)./2 + X_ellip.*sin(w_ellip) + Y_ellip.*cos(w_ellip);
% plot(x_ellip(:,:,t),y_ellip(:,:,t),'w-');


for t=1:sp:t_max
    a=a+1;
    if a == 4
        b=0; c=1;
    end
    sub_plot = subplot(2,3,a);
    set(sub_plot, 'position', [0.025+(0.325*b) 0.55-(0.5*c) 0.3 0.4]);
    pcolor(x, y, theta_rad(:,:,t)*180/pi);
    colormap Pink; shading interp;
    axis image; axis off;
    colorbar;
    set(gca,'CLim',[0 180]);
    set(gca,'Ydir','reverse');
    hold on;
    hv = quiver(xq, yq, DispXq_bd(:,:,t),DispYq_bd(:,:,t));
    set(hv,'linewidth',1, 'Color', 'k');
    hq = quiver(xq, yq, DispXq_ph(:,:,t)*pix_size,DispYq_ph(:,:,t)*pix_size);
    set(hq,'linewidth',1, 'Color', 'b');
    if t_num < 30
        title({'Vel(blk)-Tr(blu) angle distribution', ['(Degree, ' num2str((t-1)*10) 'min)']});
    elseif t_num >= 30
        title({'Vel(blk)-Tr(blu) angle distribution', ['(Degree, ' num2str(round((t-1)/6)) ' hours)']});
    end
    b = b+1;
    if a == 6
        break
    end
end

print([img_dir '/Figures/Vel_trac_angle'], '-dpng','-r300','-r0');
close;

%%
% %% Average over time
% % t_num = 72;
% fig_Overtime = figure;
% set(fig_Overtime, 'color','w','position', [500 50 600 900]);
% set(gcf,'PaperPositionMode','auto');
% time_ph = [1:t_num-1].*(1/6);
% 
% for t=1:t_num-1
% %     Disp_ph_int(:,t) = reshape(Disp_ph(:,:,t), M_ph^2, 1);
% %     do = round(imresize(reshape(domain_time{t}, M_bd, M_bd),[M_ph N_ph]));
%     
%     Disp_ph_int_dom(:,:,t) = Disp_ph(:,:,t).*round(imresize(reshape(domain_time{t}, M_bd, M_bd),[M_ph N_ph]));
%     [~,~,Disp_v{t}] = find(Disp_ph_int_dom(:,:,t));
%     Disp_ph_median(t) = median(Disp_v{t}(:,:))*pix_size/10;
%     
%     Traction_int_dom(:,t) = reshape(sqrt(tx(:,:,t).^2 + ty(:,:,t).^2), M_bd^2, 1).*domain_time{t};
%     [~,~,Traction_v{t}] = find(Traction_int_dom(:,t));
%     % Traction_int_dom(:,t) = sqrt(traction_time{t}(:,3).^2 + traction_time{t}(:,4).^2);
%     Traction_median(t) = median(Traction_v{t}(:));
%     
%     Tension_int_dom(1:size((stress_time{t}),1),t) = (stress_time{t}(:,4) + stress_time{t}(:,5))/2;
%     [~,~,Tension_v{t}] = find(Tension_int_dom(:,t));
%     Tension_median(t) = median(Tension_v{t}(:));
% end
% % yaxis = round(t_num/6);
% subplot(3,1,1);
% plot(time_ph, Disp_ph_median ,'ro'); axis([0 12 0 0.4]); % tight;
% title('Speed (\mum/min)'); xlabel('hour'); ylabel('\mum/min');
% subplot(3,1,2);
% plot(time_ph, Traction_median ,'go'); axis([0 12 0 200]); % tight;
% title('Traction (Pa)'); xlabel('hour'); ylabel('Pa');
% subplot(3,1,3);
% plot(time_ph, Tension_median ,'bo'); axis([0 12 0 500]); % tight;
% title('Tension (Pa)'); xlabel('hour'); ylabel('Pa');
% 
% print([img_dir '/Figures/All_over_time'], '-dpng','-r300','-r0');
% % print([img_dir '/Figures/All_over_time'], '-depsc');
% close;


%% Trajectory of Centroid
fg2 = figure;
set(fg2,'color','w','position', [0 0 900 850]);
set(gcf,'PaperPositionMode','auto');

for t=1:t_max
    
    s = regionprops(domain_resize_ph(:,:,t),'centroid', 'Area');
    [~,index] = max([s.Area]);
    cent = cat(1, s(index).Centroid);
    scatter(cent(1), cent(2), 50, (t-1)/6, '*');
    
    if t == 1
        cent_fst = cent;
    end
    
    hold on;
end
hold off;
c = colorbar; c.Label.String = 'Time (hr)'; c.Label.FontSize = 25;
set(gca,'CLim',[0 (t_max)/6]); colormap(fg2, cool);
title('Trajectory of centroid','FontSize',25);
xlabel('\mum','fontsize',25);
ylabel('\mum','fontsize',25);
axis equal;
set(gca,'box','on','fontsize',15);
set(gca,'Ydir','reverse');
axis ([cent_fst(1)-5 cent_fst(1)+5 cent_fst(2)-5 cent_fst(2)+5]);


print([img_dir '/Figures/CoM shifting'], '-dpng','-r300','-r0');
close;

%% Theta

fig_TractionR = figure;
set(fig_TractionR,'color','w','position', [55 20 1200 1200]);
set(gcf,'PaperPositionMode','auto');
a = 1; b = 4; c = 7;
    for t=1:24:49
    sub_plot = subplot(3,3,a);
    [Theta, ~] = cart2pol(dx_ph(:,:,t),dy_ph(:,:,t));
    Theta_line = reshape(Theta(1:M_ph,1:N_ph),M_ph^2,1);
    Theta_line(Theta_line == 0) = NaN;
    [tout, rout] = rose(Theta_line); 
    relativefreq = rout/length(Theta_line);
    pp = 0 : .01 : 0.05 * pi;
    P = polar(pp, 0.05 * ones(size(pp)));
    set(P, 'Visible', 'off')
    hold on;
    polar(tout,relativefreq)
    set(gca,'YDir','reverse');%, 'XDir','reverse');%set(gca,'View',[-90 90],'YDir','reverse');
    title({['Velocity angle distribution (', num2str(round((t-1)/6)) ' hr)']});

    clear Theta Theta_line;
    a = a+1;
    end
    %axis image; axis off;

    for t=1:24:49
    sub_plot = subplot(3,3,b);
    [Theta, ~] = cart2pol(tx(:,:,t),ty(:,:,t));
    Theta_line = reshape(Theta(1:M_ph,1:N_ph),M_ph^2,1);
    Theta_line(Theta_line == 0) = NaN;
    [tout, rout] = rose(Theta_line); 
    relativefreq = rout/length(Theta_line);
    pp = 0 : .01 : 0.1 * pi;
    P = polar(pp, 0.1 * ones(size(pp)));
    set(P, 'Visible', 'off')
    hold on;
    polar(tout,relativefreq)
    set(gca,'YDir','reverse');%, 'XDir','reverse');%set(gca,'View',[-90 90],'YDir','reverse');
    title({['Traction angle distribution (', num2str(round((t-1)/6)) ' hr)']});

    clear Theta Theta_line;
    b = b+1;
    end
    
    for t=1:24:49
    sub_plot = subplot(3,3,c);
    Theta_line = reshape(theta_rad(1:M_ph,1:N_ph),M_ph^2,1);
    Theta_line(Theta_line == 0) = NaN;
    [tout, rout] = rose(Theta_line); 
    relativefreq = rout/length(Theta_line);
    pp = 0 : .01 : 0.05 * pi;
    P = polar(pp, 0.05 * ones(size(pp)));
    set(P, 'Visible', 'off')
    hold on;
    polar(tout,relativefreq)
    set(gca,'View',[0 90]);
    title({['Vel-Tr angle distribution (', num2str(round((t-1)/6)) ' hr)']});

    clear Theta Theta_line;
    c = c+1;
    end

print([img_dir '/Figures/Angle distribution'], '-dpng','-r300','-r0');
close;

%% Kymograph
% Number of radial coords (columns)
col_min_ph = size(Vel_rad_aver,1);
col_min_bd = size(Tr_rad_aver,1);
Fntsize = 12;
% Radial coordinates (um)
r_ph = (1:col_min_ph)*grid_space_ph*pix_size; 
r_bd = (1:col_min_bd)*grid_space_bd*pix_size; 

hf1 = figure;
set(hf1,'color','w','position',[450 5 400 900]);
set(gcf,'PaperPositionMode','auto');
subplot(3,1,1);
% Kymograph of radial velocity
imagesc([0, t_ph/6],[0, col_min_ph*grid_space_ph*pix_size],Vel_rad_aver);
caxis([-0.6 0.6]); 
colormap(map); hvel = colorbar; hvel.Label.String = '(\mum/min)'; hvel.FontSize = Fntsize;
axis xy; axis tight; set(gca,'box','off');
xlabel('Time (hr)','fontsize',Fntsize);
ylabel('Radial position (\mum)','fontsize',Fntsize);
set(gca,'fontsize',Fntsize);
title('Velocity_{radial} (\mum/min)','fontsize',Fntsize);

subplot(3,1,2);
imagesc([0, t_bd/6],[0, col_min_bd*grid_space_bd*pix_size],Tr_rad_aver);
caxis([-20 20]); 
colormap(map); hvel = colorbar; hvel.Label.String = '(Pa)'; hvel.FontSize = Fntsize;
axis xy; axis tight; set(gca,'box','off');
xlabel('Time (hr)','fontsize',Fntsize);
ylabel('Radial position (\mum)','fontsize',Fntsize);
set(gca,'fontsize',Fntsize);
title('Traction_{radial} (Pa)','fontsize',Fntsize);

subplot(3,1,3);
imagesc([0, t_bd/6],[0, col_min_bd*grid_space_bd*pix_size],Tn_rad_aver);
caxis([-100 100]); 
colormap(map); hvel = colorbar; hvel.Label.String = '(Pa)'; hvel.FontSize = Fntsize;
axis xy; axis tight; set(gca,'box','off');
xlabel('Time (hr)','fontsize',Fntsize);
ylabel('Radial position (\mum)','fontsize',Fntsize);
set(gca,'fontsize',Fntsize);
title('Tension (Pa)','fontsize',Fntsize);

print([img_dir '/Figures/Radial_velocity_Kymograph'], '-dpng','-r0');
close;

%% Cell image (time)

Cell_time = figure;
set(Cell_time, 'color','w','position', [500 200 640 640]);
set(gcf,'PaperPositionMode','auto');
axis image; axis off;
set(gca,'box','on','fontsize',11);
set(gca,'Ydir','reverse');
set(gca,'color','k');

hold on;
for t = 1:t_num-1
    im_t = imread([crpdimg_dir '/phase_crpt.tif'],t);
    imagesc(im_t);colormap(Cell_time,gray);axis image;axis off;

%     hold on;
%     hv = quiver(xq, yq, DispXq_bd(:,:,t),DispYq_bd(:,:,t));
%     set(hv,'linewidth',1, 'Color', 'k');
    hq = quiver(xq, yq, DispXq_ph(:,:,t)*pix_size,DispYq_ph(:,:,t)*pix_size);
    set(hq,'linewidth',2, 'Color', 'k');
    title({'Cell image', ['(' num2str((t-1)*10) 'min)']});
    
    Cell_frame = getframe(1);
    Cell_currim = frame2im(Cell_frame);
    [Cell_A, Cell_map] = rgb2ind(Cell_currim,256);
    if t == 1;
		imwrite(Cell_A,Cell_map,[img_dir '/Figures/Cell_time.gif'],'gif','LoopCount',Inf,'DelayTime',0.2);
	else
		imwrite(Cell_A,Cell_map,[img_dir '/Figures/Cell_time.gif'],'gif','WriteMode','append','DelayTime',0.2);
	end
end
close;

%% Velocity (time)

Vel_time = figure;
set(Vel_time, 'color','w','position', [500 200 640 640]);
set(gcf,'PaperPositionMode','auto');
c = colorbar; c.Label.String = '\mum/min'; c.Label.FontSize = 15;
set(gca,'CLim',[0 1]);
axis image; axis off;
set(gca,'box','on','fontsize',11);
set(gca,'Ydir','reverse');
set(gca,'color','k');

hold on;
for t = 1:t_num-1
    pcolor(x_ph,y_ph,Vel_rho(:,:,t));
    colormap jet; shading interp;

%     hold on;
%     hv = quiver(xq, yq, DispXq_bd(:,:,t),DispYq_bd(:,:,t));
%     set(hv,'linewidth',1, 'Color', 'k');
    hq = quiver(xq, yq, DispXq_ph(:,:,t)*pix_size,DispYq_ph(:,:,t)*pix_size);
    set(hq,'linewidth',2, 'Color', 'k');
    title({'Velocity', ['(\mum/min, ' num2str((t-1)*10) 'min)']});
    
    Vel_frame = getframe(1);
    Vel_currim = frame2im(Vel_frame);
    [Vel_A, Vel_map] = rgb2ind(Vel_currim,256);
    if t == 1;
		imwrite(Vel_A,Vel_map,[img_dir '/Figures/Vel_time.gif'],'gif','LoopCount',Inf,'DelayTime',0.2);
	else
		imwrite(Vel_A,Vel_map,[img_dir '/Figures/Vel_time.gif'],'gif','WriteMode','append','DelayTime',0.2);
	end
end
close;

%% Traction (time)

Trac_time = figure;
set(Trac_time, 'color','w','position', [500 200 640 640]);
set(gcf,'PaperPositionMode','auto');
c = colorbar; c.Label.String = 'Pa'; c.Label.FontSize = 15;
set(gca,'CLim',[0 50]);
axis image; axis off;
set(gca,'box','on','fontsize',11);
set(gca,'Ydir','reverse');
set(gca,'color','k');

hold on;
for t = 1:t_num-1
    pcolor(x, y, Tr_rho(:,:,t)); 
    colormap hot; shading interp;

%     hold on;
    hv = quiver(xq, yq, DispXq_bd(:,:,t),DispYq_bd(:,:,t));
    set(hv,'linewidth',2, 'Color', 'w');
%     hq = quiver(xq, yq, DispXq_ph(:,:,t)*pix_size,DispYq_ph(:,:,t)*pix_size);
%     set(hq,'linewidth',1, 'Color', 'b');
    title({'Traction', ['(Pa, ' num2str((t-1)*10) 'min)']});
    
    Trac_frame = getframe(1);
    Trac_currim = frame2im(Trac_frame);
    [Trac_A, Trac_map] = rgb2ind(Trac_currim,256);
    if t == 1;
		imwrite(Trac_A,Trac_map,[img_dir '/Figures/Traction_time.gif'],'gif','LoopCount',Inf,'DelayTime',0.2);
	else
		imwrite(Trac_A,Trac_map,[img_dir '/Figures/Traction_time.gif'],'gif','WriteMode','append','DelayTime',0.2);
	end
end
close;

%% Tension (time)

Ten_time = figure;
set(Ten_time, 'color','w','position', [500 200 640 640]);
set(gcf,'PaperPositionMode','auto');
c = colorbar; c.Label.String = 'Pa'; c.Label.FontSize = 15;
set(gca,'CLim',[0 500]);
axis image; axis off;
set(gca,'box','on','fontsize',11);
set(gca,'Ydir','reverse');
set(gca,'color','k');

hold on;
for t = 1:t_num-1
    pcolor(x, y, Tn_rho(:,:,t));
    colormap hot; shading interp;

%     hold on;
%     hv = quiver(xq, yq, DispXq_bd(:,:,t),DispYq_bd(:,:,t));
%     set(hv,'linewidth',1, 'Color', 'k');
%     hq = quiver(xq, yq, DispXq_ph(:,:,t)*pix_size,DispYq_ph(:,:,t)*pix_size);
%     set(hq,'linewidth',1, 'Color', 'b');
    title({'Tensional stress', ['(Pa, ' num2str((t-1)*10) 'min)']});
    
    Ten_frame = getframe(1);
    Ten_currim = frame2im(Ten_frame);
    [A, map] = rgb2ind(Ten_currim,256);
    if t == 1;
		imwrite(A,map,[img_dir '/Figures/Tension_time.gif'],'gif','LoopCount',Inf,'DelayTime',0.2);
	else
		imwrite(A,map,[img_dir '/Figures/Tension_time.gif'],'gif','WriteMode','append','DelayTime',0.2);
	end
end
close;
%% Average over time with histogram
t_ph = t_max;

for t = 1 : t_ph
    Rho_vel_line(:,t) = reshape(Vel_rho(1:M_ph,1:N_ph,t),M_ph^2,1);
    Rho_tr_line(:,t) = reshape(Tr_rho(1:M_bd,1:N_bd,t),M_bd^2,1);
    Rho_tn_line(:,t) = reshape(Tn_rho(1:M_bd,1:N_bd,t),M_bd^2,1);
    Rho_sh_line(:,t) = reshape(sh(1:M_bd,1:N_bd,t),M_bd^2,1);
    Rho_rs_line(:,t) = reshape(Tn_rho(1:M_bd,1:N_bd,t)./sh(1:M_bd,1:N_bd,t),M_bd^2,1);
    Rho_th_line(:,t) = reshape(theta_ellips(1:M_ph,1:N_ph,t),M_ph^2,1);
end

Rho_vel_line(Rho_vel_line == 0) = NaN;
Rho_tr_line(Rho_tr_line == 0) = NaN;
Rho_tn_line(Rho_tn_line == 0) = NaN;
Rho_sh_line(Rho_sh_line == 0) = NaN;
Rho_rs_line(Rho_rs_line == 0) = NaN;
Rho_th_line(Rho_th_line == 0) = NaN;
%Rho_th_line = Rho_th_line.*(pi/180);

Edge_min_vel = floor(min(min(Rho_vel_line))*0.1)*10;
Edge_max_vel = 1; %0.8; %ceil(max(max(Rho_ph_line))*0.1)*10;
Edge_vel = linspace(Edge_min_vel, Edge_max_vel, 50);

Edge_min_tr = floor(min(min(Rho_tr_line))*0.1)*10;
Edge_max_tr = 40; %ceil(max(max(Rho_ph_line))*0.1)*10;
Edge_tr = linspace(Edge_min_tr, Edge_max_tr, 50);

Edge_min_tn = floor(min(min(Rho_tn_line))*0.1)*10;
Edge_max_tn = 600; %ceil(max(max(Rho_ph_line))*0.1)*10;
Edge_tn = linspace(Edge_min_tn, Edge_max_tn, 50);

Edge_min_sh = floor(min(min(Rho_sh_line))*0.1)*10;
Edge_max_sh = 200; %ceil(max(max(Rho_ph_line))*0.1)*10;
Edge_sh = linspace(Edge_min_sh, Edge_max_sh, 50);

Edge_min_rs = floor(min(min(Rho_rs_line))*0.1)*10;
Edge_max_rs = 10; %ceil(max(max(Rho_ph_line))*0.1)*10;
Edge_rs = linspace(Edge_min_rs, Edge_max_rs, 50);


Edge_min_th = floor(min(min(Rho_th_line))*0.1)*10;
Edge_max_th = 180; %ceil(max(max(Rho_ph_line))*0.1)*10;
Edge_th = linspace(Edge_min_th, Edge_max_th, 50);

for t = 1 : t_ph
    Num_vel(:,t) = histcounts(Rho_vel_line(:,t), Edge_vel, 'Normalization', 'probability');
    Num_tr(:,t) = histcounts(Rho_tr_line(:,t), Edge_tr, 'Normalization', 'probability');
    Num_tn(:,t) = histcounts(Rho_tn_line(:,t), Edge_tn, 'Normalization', 'probability');
    Num_sh(:,t) = histcounts(Rho_sh_line(:,t), Edge_sh, 'Normalization', 'probability');
    Num_rs(:,t) = histcounts(Rho_rs_line(:,t), Edge_rs, 'Normalization', 'probability');
    Num_th(:,t) = histcounts(Rho_th_line(:,t), Edge_th, 'Normalization', 'probability');
end

fig_HistoAverage = figure;
set(fig_HistoAverage, 'color','w','position', [500 5 1200 1200]);
set(gcf,'PaperPositionMode','auto');
time_ph = [1:t_num-1].*(1/6);
Fntsize = 15;

av1 = subplot(3,2,1);
t_all = linspace(0, t_ph/6, t_ph);
% Kymograph of radial velocity
imagesc(t_all,[0, Edge_max_vel],Num_vel*100);
colormap(av1, hot);
caxis([0 15]); 
colormap hot; hvel = colorbar; % cbfreeze(hc);
hvel.Label.String = 'Probability (%)';
hvel.FontSize = Fntsize;
axis xy; axis ([0 t_ph/6 0 Edge_max_vel]); set(gca,'box','off');
xlabel('Time (hr)','fontsize',Fntsize);
ylabel('Average speed (\mum/min)','fontsize',Fntsize);
set(gca,'fontsize',Fntsize);
title('Average speed_{histogram} (\mum/min)','fontsize',Fntsize);
hold on
Rho_vel_aver = median(Rho_vel_line, 'omitnan');
plot(t_all, Rho_vel_aver ,'wo');

av3 = subplot(3,2,3);
t_all = linspace(0, t_ph/6, t_ph);
% Kymograph of radial velocity
imagesc(t_all,[0, Edge_max_tr],Num_tr*100);
colormap(av3, hot);
caxis([0 10]); 
colormap cool; htr = colorbar; % cbfreeze(hc);
htr.Label.String = 'Probability (%)';
htr.FontSize = Fntsize;
axis xy; axis ([0 t_ph/6 0 Edge_max_tr]); set(gca,'box','off');
xlabel('Time (hr)','fontsize',Fntsize);
ylabel('Average traction (Pa)','fontsize',Fntsize);
set(gca,'fontsize',Fntsize);
title('Average traction_{histogram} (Pa)','fontsize',Fntsize);
hold on
Rho_tr_aver = median(Rho_tr_line, 'omitnan');
plot(t_all, Rho_tr_aver,'wo');

av5 = subplot(3,2,5);
t_all = linspace(0, t_ph/6, t_ph);
% Kymograph of radial velocity
imagesc(t_all,[0, Edge_max_tn],Num_tn*100);
colormap(av5, hot);
caxis([0 15]); 
colormap hot; htn = colorbar; % cbfreeze(hc);
htn.Label.String = 'Probability (%)';
htn.FontSize = Fntsize;
axis xy; axis ([0 t_ph/6 0 Edge_max_tn]); set(gca,'box','off');
xlabel('Time (hr)','fontsize',Fntsize);
ylabel('Average tension (Pa)','fontsize',Fntsize);
set(gca,'fontsize',Fntsize);
title('Average tension_{histogram} (Pa)','fontsize',Fntsize);
hold on
Rho_tn_aver = median(Rho_tn_line, 'omitnan');
plot(t_all, Rho_tn_aver ,'wo');

av4 = subplot(3,2,4);
t_all = linspace(0, t_ph/6, t_ph);
% Kymograph of radial velocity
imagesc(t_all,[0, Edge_max_sh],Num_sh*100);
colormap(av4, hot);
caxis([0 10]); 
colormap hot; htn = colorbar; % cbfreeze(hc);
htn.Label.String = 'Probability (%)';
htn.FontSize = Fntsize;
axis xy; axis ([0 t_ph/6 0 Edge_max_sh]); set(gca,'box','off');
xlabel('Time (hr)','fontsize',Fntsize);
ylabel('Max shear stress (Pa)','fontsize',Fntsize);
set(gca,'fontsize',Fntsize);
title('Max shear stress_{histogram} (Pa)','fontsize',Fntsize);
hold on
Rho_sh_aver = median(Rho_sh_line, 'omitnan');
plot(t_all, Rho_sh_aver ,'wo');

av6 = subplot(3,2,6);
t_all = linspace(0, t_ph/6, t_ph);
% Kymograph of radial velocity
imagesc(t_all,[0, Edge_max_rs],Num_rs*100);
colormap(av6, hot);
caxis([0 15]); 
colormap hot; htn = colorbar; % cbfreeze(hc);
htn.Label.String = 'Probability (%)';
htn.FontSize = Fntsize;
axis xy; axis ([0 t_ph/6 0 Edge_max_rs]); set(gca,'box','off');
xlabel('Time (hr)','fontsize',Fntsize);
ylabel('Ratio','fontsize',Fntsize);
set(gca,'fontsize',Fntsize);
title('Tension/Max shear stress_{histogram} (R)','fontsize',Fntsize);
hold on
Rho_rs_aver = median(Rho_rs_line, 'omitnan');
plot(t_all, Rho_rs_aver ,'wo');


av2 = subplot(3,2,2);
t_all = linspace(0, t_ph/6, t_ph);
% Kymograph of radial velocity
imagesc(t_all,[0, 180],Num_th*100);
colormap(av2, hot);
caxis([0 8]); 
colormap hot; htn = colorbar; % cbfreeze(hc);
htn.Label.String = 'Probability (%)';
htn.FontSize = Fntsize;
axis xy; axis ([0 t_ph/6 0 180]); set(gca,'box','off');
xlabel('Time (hr)','fontsize',Fntsize);
ylabel('Degree','fontsize',Fntsize);
set(gca,'fontsize',Fntsize);
title('Theta(Vel-Tr)_{histogram} (Degree)','fontsize',Fntsize);
hold on
Rho_th_aver = median(Rho_th_line, 'omitnan');
plot(t_all, Rho_th_aver ,'wo');

print([img_dir '/Figures/Average_all_histogram'], '-dpng','-r300','-r0');
close;
end



