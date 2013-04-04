% Tutorial DCM

clear; close all; clc;
cfg = [];
addpath(genpath( 'D:/libraries/matlab/fieldtrip'))
ft_defaults();

subject  = 'dcm_data/subject1';

cfg.dataset     =   strcat(subject,'.bdf');
cfg.continuous  =   'yes';
data_org        =   ft_preprocessing(cfg);

clear cfg;
cfg.channel =       'all';
cfg.reref =         'yes';
cfg.refchannel =    'all';
% cfg.refchannel =    {'all'};
cfg.lpfilter =      'yes';
cfg.hpfilter =      'yes';
cfg.lpfreq   =      30;
cfg.hpfreq   =      0.5;
cfg.demean   =      'yes';
data_org        = ft_preprocessing(cfg,data_org);