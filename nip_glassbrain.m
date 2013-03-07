function h = nip_glassbrain(dip_pos, J, a_handle)
% h = nip_glassbrain(dip_pos, J, a_handle)
% Shows the activity 'J' of the dipoles in the dipoles position dip_pos
% Input:
%       dip_pos -> Ndx3. 3D coordinates of the position of the dipoles
%       J -> Ndx1. Vector with the activity on the dipoles
%       a_handle -> axes handle. This is where the plot is going to be
%               drawn.
% Output:
%       h -> axes handle.
%
% Additional comments :
% 1) This function is a Wrapper to the spm_mip function of the SPM8
% software
% 2) Do not use this function if the data used (cortex mesh, lead field etc...) is the one provided in this
% toolbox, to use it, you should rotate it in order to fit the template
% provided in the spm toolbox). Use nip_reconstruction3d  instead
% 
% Juan S. Castano C.
% jscastanoc@gmail.com
% 15 Feb 2013

axes(a_handle);
min_power = 2e-1;
E = J.^2;
E = E - min(min(E));
E = E/max(E);
j= find(E > min_power);
E = E(j);
voxel_sz = 10;
spm_mip(E,dip_pos(j,:),voxel_sz); 
colormap gray
axis image
h = gca;
end
