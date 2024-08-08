function stress

clear;
close all;
clc;

% Load data
load('traction.mat');
load('inp_setup.mat');
% 
% % Total of correlations
% t_num = size(tx,3);
% Make directories for FE based stress computation code
mkdir('Traction'); mkdir('FEM');

% Preallocate cells to save nodes and stresses when running FE
node_coords_time = cell(t_num,1);
stress_time = node_coords_time;
domain_time = node_coords_time;
traction_time = node_coords_time;
[M, N] = size(x);

for t = 1:t_num
    % Get t-th tractions
    tx_t = tx(:,:,t);
    ty_t = ty(:,:,t);
    % Get t-th domain
    domain = imread([crpdimg_dir '/domain.tif'],t);
%    domain = double(domain);
    domain = imresize(domain, [M N]);
%    domain = round(domain);
    if max(domain(:))>0 % Only run analysis if there are nonzero values of domain (ie, the island was found by the edge finding algorithm)
%         domain = domain/max(domain(:)); % Normalize to 1
%         
%         % Smooth domain again
        SE = strel('disk',2,0);
        domain = imclose(domain,SE); % Morphological closing is same as dilation and then erosion
        domain = imdilate(domain,SE);
        s = regionprops(domain,'Area');
        domain = bwareafilt(domain,[max([s.Area])-1 max([s.Area])]);
        domain = imfill(domain,'holes');
        domain = imclose(domain,SE);
%         
        % --- TRACTION/MOMENT CORRECTION ---
        % Equation for computing moment induced by tractions about center of island
%         moment = @(t_shift) sum(sum( (x-xc).*(ty_domain-t_shift) - (y-yc).*(tx_domain-t_shift) ));
%         % Compute shift required in traction data by setting moment equal to zero
%         t_shift = fzero(moment,0);
        % Moment corrected traction data. (Note, this is still not in force
        % balance!)
        tx_mzero = tx_t ;%- t_shift;
        ty_mzero = ty_t ;%- t_shift;
        % Correct for force balance
        tx_mzero(domain==1) = tx_mzero(domain==1) - mean(tx_mzero(domain==1));
        ty_mzero(domain==1) = ty_mzero(domain==1) - mean(ty_mzero(domain==1));
        % Set tractions outside domain to zero (this isn't strictly
        % necessary for stress computation, but it makes the size of the 
        % dat files smaller)
        tx_mzero(domain==0) = 0;
        ty_mzero(domain==0) = 0;
        % --------------------------------------------
        
        % Assemble traction array to save
        %     traction_array = [x(:) y(:) tx_t(:) ty_t(:)];
        traction_array = [x(:) y(:) tx_mzero(:) ty_mzero(:)];
        % Convert domain to vector to save
        domain_vector = double(domain(:));
        
        % Delete files from previous loops
        delete('FEM/bforce.dat'); delete('FEM/node_coords.dat'); delete('FEM/StressStrain.dat');
        % Save data
        save('Traction/domain.dat','domain_vector','-ascii');
        save('Traction/traction.dat','traction_array','-ascii');
        
        % Run FE based stress computation code
        % Note this can be called using unix('island.exe') or
        % system('island.exe'), but using the bang command allows for using the
        % & symbol to tell Matlab to return to the script and continue running
        % without requiring user to hit enter to continue.
        
        !I:\Analysing/island.exe&
        
        % Since I use the & symbol, I need to wait for Windows to finish
        % running
        pause(3); % Pause for a few seconds
        % system([code_dir '/island.exe'])
        
        % Load compuated stress data
        domain_time{t} = load('Traction/domain.dat'); % Free edges of cell monolayer
        % Format: column of binary (1 or 0)
        traction_time{t} = load('Traction/traction.dat');
        % Format: x | y | tx | ty
        
        if exist('FEM/StressStrain.dat','file') == 2 % exist will equal 2 if looking for a file
            node_coords_t = load('FEM/node_coords.dat');
            stress_t = load('FEM/StressStrain.dat');
            node_coords_time{t} = node_coords_t; % Need to store as a cell because domain is different for each time point and therefore number of nodes is different for each time point
            stress_time{t} = stress_t;
            disp(['Stress computation for time ',num2str(t),' of ',num2str(t_num),' complete.'])
        else
            disp(['Unable to compute stresses for time ',num2str(t)'.'])
        end
        
    else
        disp(['Unable to compute stresses for time ',num2str(t)'.'])
    end
    
end

% Save data
save('stresses.mat','domain_time','traction_time','node_coords_time','stress_time');

% Remove directories for FE based stress computation code
rmdir('Traction','s'); rmdir('FEM','s');
end

