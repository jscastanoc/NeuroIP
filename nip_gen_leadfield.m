function lf = nip_gen_leadfield(mesh_head, dip_pos, elec, vol)
% lf = nip_gen_leadfield(mesh_head, dip_pos, elec, vol)
% Input ->
%           mesh_head   -> struct. containing the meshs that describe the
%                       geometry of the head (skin, skull, brain) or (skin);
%           dip_pos     -> matrix. Ndx3 containg the dipole positions
%           elec        -> struct. Electrode description according to the Fieldtrip
%                       toolbox:    elec.label, elec.chanpos elec.elecpos
%                       ref: http://fieldtrip.fcdonders.nl/faq/how_are_electrodes_magnetometers_or_gradiometers_described
%           vol         -> Char.
%                       .mat file saved with the info generated by vol = ft_headmodel_dipoli(head)
%                       or similar (see the fieldtrip doc for further information
%                       In the mean time, you can just use the data provided under the
%                       dir  "extra_data".
%                       if .mat file does not exist, it computes the volume
%                       conductor model with the information provided.
% output ->
%           lf          -> matrx. Ncx3Nd Lead field matrix
%
% Additional Comments:
% This function is a wrapper for the Fieldtrip funcions:
%    - ft_headmodel_dipoli - ft_prepare_vol_sens - ft_compute_leadfield.
% Juan S. Castaño C. 
% 15 Jan 2013

home = fileparts(which('nip_gen_leadfield'));
addpath(strcat(home,'/external/fieldtrip/forward'));

vol_aux = mesh_head;
head = [];
for i = 1:numel(vol_aux)
   head(i).pnt = vol_aux{i}.vertices; 
   head(i).tri = vol_aux{i}.faces;
end

if exist(fullfile(home,vol),'file')
    load(fullfile(home,vol))    
    [vol, sens] = ft_prepare_vol_sens(vol,elec);
else
    disp('Computing Volume Conductor Model')
    file_name = fullfile(home,vol);
    try
        vol = ft_headmodel_dipoli(head);
    catch
%         vol = ft_headmodel_openmeeg(head);
        vol = ft_headmodel_bemcp(head);
    end
    save(file_name,'vol');
end

[vol, sens] = ft_prepare_vol_sens(vol,elec);
disp('Computing lead fields')
lf =  ft_compute_leadfield(dip_pos, sens, vol);


