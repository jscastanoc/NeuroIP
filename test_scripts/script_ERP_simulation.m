% Simulation ERP (wavelets and biological noise)
close all; clc; clear

addpath('../external/source_toolbox/haufe/')
addpath('../external/source_toolbox/nolte/')
addpath('../external/source_toolbox/simulations/')
rng('default')
rng('shuffle');

warning off

load clab_example;
load clab_10_10;
clab = clab_10_10;
% sa = prepare_sourceanalysis(clab, 'icbm152b_sym');
sa = prepare_sourceanalysis(clab, 'montreal');


for i = 1:3
    temp(:,:) = sa.V_cortex(:,:,1);
end
L = reshape(permute(sa.V_cortex,[1 3 2]), size(sa.V_cortex,1), size(sa.V_cortex,2)*3);

cfg.fs = 200; % Sample frequency (Hz)
cfg.t = 0:1/cfg.fs:2; % Time vector (seconds)
t = cfg.t;
cfg.L = L;

model = nip_create_model(cfg);


phase_shift = [0.5 1.5] ; % Phase shift for the sources (in seconds)
Nact = length(phase_shift);
fc_wl = 5; % Central frequency for the wavelet
for i = 1:Nact
    f0 = fc_wl;
    
    % Normalization terms for the wavelet
    sigma_f = f0/7;
    sigma_t = 1/(2*pi*sigma_f);
    
    % "source" contains the time courses of active sources
    source(i,:) =  real(exp(2*1i*pi*f0*t).*...
        exp((-(t-phase_shift(i)).^2)/(2*sigma_t^2)));
end


[y, Jclean] = nip_simtrials(model.L,sa.cortex, ...
        source, model.t, 500 , 200, 10, -2);
plot(model.t,y)

