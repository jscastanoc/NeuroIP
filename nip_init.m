function status = nip_init(dir)
% status = nip_init()
% Initializes paths of the toolbox
%
% Input: dir->location where the external libraries are installed
%
% Additional Comments:
% Juan S. Castaño C.
% 19 Feb 2013

addpath(strcat(dir,'/toolbox_graph'));

addpath(strcat(dir,'/dal-master'));
addpath(strcat(dir,'/pymex-master'));
addpath(strcat(dir,'/ltfat'));

 dir = strcat(fileparts(which('nip_init')));
addpath(strcat(dir,'/help_scripts/'));
addpath(strcat(dir,'/data'));

ltfatstart;
end
