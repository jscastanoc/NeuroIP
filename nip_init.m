function status = nip_init()
% status = nip_init()
% This function is under construction, currently, it only adds some
% directories to the path, just call it before doing anything.
%
% Additional Comments:
% Juan S. Castaño C.
% 19 Feb 2013
dir = strcat(fileparts(which('nip_init')));
addpath(strcat(dir,'/data'));
addpath(strcat(dir,'/external/toolbox_graph'));
addpath(strcat(dir,'/external/toolbox_graph/toolbox'));
addpath(strcat(dir,'/external/fieldtrip'));
addpath(strcat(dir,'/external/fieldtrip/forward'));
addpath(strcat(dir,'/external/spm8'));
addpath(strcat(dir,'/external/nway310/ver3.1'));
addpath(strcat(dir,'/external/dal_ver1.05'));
addpath(strcat(dir,'/external/matlab_bgl'));
addpath(strcat(dir,'/external/FastEMD-3'));
addpath(strcat(dir,'/external/source_toolbox/haufe/'));
addpath(strcat(dir,'/external/source_toolbox/nolte/'));
addpath(strcat(dir,'/external/source_toolbox/simulations/'));
addpath(strcat(dir,'/external/bbci/toolbox/startup/'));
addpath(strcat(dir,'/external/ltfat/'));
addpath(strcat(dir,'/help_scripts/'));
addpath(strcat(dir,'/external/L1_homotopy/'));
addpath(strcat(dir,'/external/RWL1DF/'));
addpath(strcat(dir,'/external/kkmeans/'));
addpath(strcat(dir,'/external/pymex/'));

dir = '/mnt/data/Master_Results';
addpath(strcat(dir,'/external/toolbox_graph'));
addpath(strcat(dir,'/external/toolbox_graph/toolbox'));
addpath(strcat(dir,'/external/fieldtrip'));
addpath(strcat(dir,'/external/fieldtrip/forward'));
addpath(strcat(dir,'/external/spm8'));
addpath(strcat(dir,'/external/nway310/ver3.1'));
addpath(strcat(dir,'/external/dal_ver1.05'));
addpath(strcat(dir,'/external/matlab_bgl'));
addpath(strcat(dir,'/external/FastEMD-3'));
addpath(strcat(dir,'/external/source_toolbox/haufe/'));
addpath(strcat(dir,'/external/source_toolbox/nolte/'));
addpath(strcat(dir,'/external/source_toolbox/simulations/'));
addpath(strcat(dir,'/external/bbci/toolbox/startup/'));
addpath(strcat(dir,'/external/ltfat/'));
addpath(strcat(dir,'/help_scripts/'));
addpath(strcat(dir,'/external/L1_homotopy/'));
addpath(strcat(dir,'/external/RWL1DF/'));
addpath(strcat(dir,'/external/kkmeans/'));
addpath(strcat(dir,'/external/pymex/'));
addpath(strcat(dir,'/external/matgrid/'));
% ltfatstart;
% setup_path;
end
