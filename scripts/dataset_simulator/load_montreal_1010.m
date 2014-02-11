% Load Lead Field Montreal 10-10
load sa_montreal
load clab_example
load clab_10_10; % Labels for 59 channels under 10 10 protocol
clab = clab_10_10;

temp = sa.V_cortex_coarse;
L = nip_translf(temp); % Leadfield matrix
L = L(find(ismember(clab_example,clab)),:);
clear temp