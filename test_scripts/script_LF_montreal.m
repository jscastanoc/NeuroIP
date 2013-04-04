% Script to generate the lead field using the data provided in:
% http://eeg.pl/Members/jarekz/lead-field-data-for-subjects-form-the-paper/view
clear; close all; clc;
load data/sa-montreal;

full_labels = sa.clab_electrodes;


% Labels for the 10-10 standard (34 electrodes)
eeg_std = '10-10';
labels={'Cz','Oz','PO3','PO4','O1','O2','P3','P4','P8', ...
           'P7','C3','T7','C1','C2','C4','T8','Fz','F4','F8','F3','F7', ...
           'AF3','AF4','Fp1','FC5','FC1','FC2','FC6','CP5','CP1','CP2', ...
           'CP6','Fp2','Pz'};

% Labels for the emotiv (14 electrodes)
% eeg_std = 'emotiv'
% labels = {'AF3', 'F7','F3', 'FC5', 'T7', 'P7', 'O1', 'O2',...
%     'P8', 'T8', 'FC6', 'F4', 'F8', 'AF4'};


[~, idx_label] = ismember(labels,full_labels);


% cortex_mesh describing the geometry of the cortex (select the number of dipoles to consider)
cortex_mesh.vertices = sa.cortex.vc;
cortex_mesh.faces = sa.cortex.tri;
cortex_mesh = nip_subsample_mesh(cortex_mesh,1999);


% cortex_mesh describing the geometry of the head
aux_vol = sa.vc(3); % Number of shells to use 
                    % You can use the three of them but if you will need a
                    % lot of RAM and a lot of time.
for i=1:numel(aux_vol)
    head{i}.vertices = aux_vol{i}.vc;
    head{i}.faces = aux_vol{i}.tri;
end


% Sensor description for the Fieldtrip Toolbox.

elec.label = labels;
elec.chanpos = sa.locs_3D(idx_label,1:3);
elec.elecpos = sa.locs_3D(idx_label,1:3);


lf = nip_gen_leadfield(head, cortex_mesh.vertices, elec, 'extra_data/vol_dipoli_1shell');


% We assume that the orientation of the dipoles is perpendicular to the
% cortex
cortex_mesh_normals = spm_eeg_inv_normals(cortex_mesh.vertices,cortex_mesh.faces);
L = zeros(size(lf, 1), size(lf, 2)/3);
for i = 1:size(L,2)
    L(:, i) = lf(:, (3*i- 2):(3*i))*cortex_mesh_normals(i, :)';
end

clear elec_pos full_labels i idx_label labels lf cortex_mesh_normals sa vol aux_vol


save(strcat('data/montreal',num2str(size(L,2)),'_',eeg_std))