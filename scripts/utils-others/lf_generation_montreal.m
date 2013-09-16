% Script to generate the lead field using the data provided in:
% http://eeg.pl/Members/jarekz/lead-field-data-for-subjects-form-the-paper/view
clear; close all; clc;
nip_init();

% Load data
dir = '/mnt/data/Datasets/';
load(strcat(dir,'sa-montreal'));

full_labels = sa.clab_electrodes;

% Labels for the 10-10 standard (34 electrodes)
eeg_std = '10-10';
labels={'Cz','Oz','PO3','PO4','O1','O2','P3','P4','P8', ...
           'P7','C3','T7','C1','C2','C4','T8','Fz','F4','F8','F3','F7', ...
           'AF3','AF4','Fp1','FC5','FC1','FC2','FC6','CP5','CP1','CP2', ...
           'CP6','Fp2','Pz'};

% Labels for the emotiv (14 electrodes)
eeg_std = 'emotiv'
labels = {'AF3', 'F7','F3', 'FC5', 'T7', 'P7', 'O1', 'O2',...
    'P8', 'T8', 'FC6', 'F4', 'F8', 'AF4'};


[~, idx_label] = ismember(labels,full_labels);


% cortex_mesh describing the geometry of the cortex (select the number of dipoles)
cortex_mesh.vertices = sa.cortex.vc;
cortex_mesh.faces = sa.cortex.tri;
num_dip = 3999;
cortex_mesh = nip_subsample_mesh(cortex_mesh,num_dip);
Nd = size(cortex_mesh.vertices,1);

% cortex_mesh describing the geometry of the head
aux_vol = sa.vc(1:3); % Number of shells to use 
                    % You can use the three of them but if you will need a
                    % lot of RAM and a lot of time.
graphic_mode = false;
if graphic_mode        
    fig_vol = figure('Units','normalized','position',[0.2 0.2 0.5 0.5]);
end

for i=numel(aux_vol):-1:1
    head{i}.vertices = aux_vol{i}.vc;
    head{i}.faces = aux_vol{i}.tri;
    if graphic_mode
        figure(fig_vol)
        h = patch('Faces', head{i}.faces, 'Vertices', head{i}.vertices);
        axis equal
        axis off
        set(h,'facealpha',0.1);
        pause(1)
        hold on
    end
end
% break
%%
% Load precomputed model if you already have it, because this takes a lot
% of time, you have a precomputed model for this dataset under data/fp/montreal_vol.mat

vol_aux = head;
head_dip = [];
for i = 1:numel(vol_aux)
   head_dip(i).pnt = vol_aux{i}.vertices; 
   head_dip(i).tri = vol_aux{i}.faces;
end
clear head vol_aux
% vol = ft_headmodel_dipoli(head_dip);
vol = ft_headmodel_bemcp(head_dip);


%%
% Sensor description for the Fieldtrip Toolbox.
elec.label = sa.clab_electrodes;
elec.chanpos = sa.locs_3D(:,1:3);
elec.elecpos = sa.locs_3D(:,1:3);
  
[vol, sens] = ft_prepare_vol_sens(vol,elec);

L =  ft_compute_leadfield(cortex_mesh.vertices, sens, vol);

% If you want to assume that the orientation of the dipoles is perpendicular to the
% cortex, uncomment the following
% cortex_mesh_normals = spm_eeg_inv_normals(cortex_mesh.vertices,cortex_mesh.faces);
% L = zeros(size(lf, 1), size(lf, 2)/3);
% for i = 1:size(L,2)
%     L(:, i) = lf(:, (3*i- 2):(3*i))*cortex_mesh_normals(i, :)';
% end
% Lfull = lf;
% clear elec_pos full_labels i idx_label labels lf cortex_mesh_normals sa vol aux_vol Lfull


% save(fullfile(home,strcat('data/montreal',num2str(size(L,2)),'_',eeg_std)))