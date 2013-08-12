% Plot cortex and head surface with sensor labels.
clear; close all; clc;
nip_init();

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preproceso / simulacion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero de dipolos a considerar
Nd = 8000; 

% Cargar datos(lead field, mesh del cerebro etc...
load(strcat('../data/const_or/montreal',num2str(Nd),'_10-10.mat'))

vol = cortex_mesh;

data = zeros(Nd*3,1);
data(1) = 1;

h=figure;
% nip_reconstruction3d(vol,data,gca);
vol.vc = vol.vertices;
vol.tri = vol.faces;
showsurface(vol)
hold on
scatter3(elec.chanpos(:,1),...
        elec.chanpos(:,2),...
        elec.chanpos(:,3),'.')
for i=1:size(L,1)
    if strcmp(elec.label{i},'Fp2') || strcmp(elec.label{i},'Fp1')
        yoffset = -1;
    else
        yoffset = 0;
    end
    text(elec.chanpos(i,1)+0.5,...
        elec.chanpos(i,2)+yoffset,...
        elec.chanpos(i,3),elec.label{i},'FontSize',14)
end
%savefig('sensors',h,'eps')
