% Script to reconstruct EEG from the eeg.pl dataset
nip_init();
%% Load EEG and select window of interest
clear; close all; clc;
cfg = [];

subject  = 'FRAANN';

cfg.dataset     =   strcat(subject,'_EEG_DATA.edf');
cfg.continuous  =   'yes';
data_org        =   ft_preprocessing(cfg);

clear cfg;
cfg.channel =       'all';
cfg.reref =         'yes';
cfg.refchannel =    'all';
cfg.lpfilter =      'yes';
cfg.hpfilter =      'yes';
cfg.lpfreq   =      25;
cfg.hpfreq   =      0.5;
% cfg.demean   =      'yes';
data_org        = ft_preprocessing(cfg,data_org);

% time_window = data_org.fsample*[6*60+56 6*60+57];
time_window = data_org.fsample*[45 46];

t = data_org.time{1}(time_window(1):time_window(2));
y = data_org.trial{1}(:, time_window(1):time_window(2));


% plot(t, y(1,:))
% xlabel('time (s)')
% ylabel('channel amplitude (a.u.)')
% title(sprintf('trial %d', trialsel));


%% Load Lead field 

% load leadfield;
load FRAANN.mat;
clear cfg;
% cfg.L = lf/max(abs(lf(:)));
cfg.L = lf;
clear lf
cfg.fs = data_org.fsample; % Frecuencia de muestreo para la simulacion
cfg.t = t; % Vector de tiempo
cfg.y = y;
cfg.cortex = src;
clear y t;
model = nip_create_model(cfg);


labels = elec.label;


% There are different standards, the ones involved in the current test
% change the name of the electrodes like this:
% four electrodes have different names compared to the 10-20 system; 
% these are T7, T8, P7, and P8. 
% T6 -> P8
% T4 -> T8
% T3 -> T7
% T5 -> P7
% ref: http://www.bem.fi/book/13/13.htm#03
label_change = {{'P7','T5'},{'T7','T3'},{'T8','T4'},{'P8','T6'}};

for i = 1:numel(label_change)
    cr_pair = label_change{i};
    ch = ismember(labels,cr_pair{1});
    ch = find(ch);
    labels{ch} = cr_pair{2};
end
sel_ch = find(ismember(labels,data_org.label));

% Reduce lead field matrix according to the available channels
model.y = model.y(find(ismember(labels(sel_ch),data_org.label)),:);
model.L = model.L(sel_ch,:);


Q = diag(nip_lcmv(model.y,model.L));
Q = Q/max(Q(:));

model.y = model.y/max(abs(model.y(:)));

[J_rec,~] = nip_loreta(model.y,model.L,Q);

figure;
nip_reconstruction3d(model.cortex, diag(Q) ,struct('axes',gca))