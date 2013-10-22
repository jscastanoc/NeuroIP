function status = nip_init()
% status = nip_init()
% This function is under construction, currently, it only adds some
% directories to the path, just call it before doing anything.
%
% Additional Comments:
% Juan S. Castaño C.
% 19 Feb 2013
addpath(strcat(fileparts(which('nip_init'))));
addpath(strcat(fileparts(which('nip_init')),'/data'));
addpath(strcat(fileparts(which('nip_init')),'/external/toolbox_graph'));
addpath(strcat(fileparts(which('nip_init')),'/external/toolbox_graph/toolbox'));
addpath(strcat(fileparts(which('nip_init')),'/external/fieldtrip'));
addpath(strcat(fileparts(which('nip_init')),'/external/fieldtrip/forward'));
addpath(strcat(fileparts(which('nip_init')),'/external/spm8'));
addpath(strcat(fileparts(which('nip_init')),'/external/nway310/ver3.1'));
addpath(strcat(fileparts(which('nip_init')),'/external/dal_ver1.05'));
addpath(strcat(fileparts(which('nip_init')),'/external/matlab_bgl'));
addpath(strcat(fileparts(which('nip_init')),'/external/FastEMD-3'));
addpath(strcat(fileparts(which('nip_init')),'/external/source_toolbox/haufe/'));
addpath(strcat(fileparts(which('nip_init')),'/external/source_toolbox/nolte/'));
addpath(strcat(fileparts(which('nip_init')),'/external/source_toolbox/simulations/'));
addpath(strcat(fileparts(which('nip_init')),'/external/bbci/toolbox/startup/'));
addpath(strcat(fileparts(which('nip_init')),'/external/ltfat/'));
addpath(strcat(fileparts(which('nip_init')),'/help_scripts/'));
addpath(strcat(fileparts(which('nip_init')),'/external/L1_homotopy/'));
addpath(strcat(fileparts(which('nip_init')),'/external/RWL1DF/'));
addpath(strcat(fileparts(which('nip_init')),'/external/kkmeans/'));
ltfatstart;
setup_path;
end
