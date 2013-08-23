% session_list = get_session_list('projects/AudioVisualSpeller');
clc; clear; close;
startup_bbci;

session_list = {'VPnal_12_11_15','VPnam_12_11_15','VPnan_12_11_16','VPnao_12_11_19',...
'VPnap_12_11_20','VPnaq_12_11_22','VPnar_12_11_23','VPnas_12_11_28','VPnat_12_11_29',...
'VPnau_12_11_30','VPhbo_12_12_01','VPlg_12_12_05','VPhbn_12_12_06','VPnav_12_12_08',...
'VPnaw_12_12_11'}

stim_target = [31:36 81:86 111:116];
stim_nontarget = [11:16 61:66 91:96];

ii=13 %for ii = 1:length(session_list)
icond = 1;
% dir = 'D:/Datasets/real_data/';
fname = strcat(session_list{ii}, '/AudiVisual_Depend_', num2str(icond) , '*');



[epo, epo_r] = stdERPanalysis(fname, {stim_target stim_nontarget; 'Target' 'Non-Target'}, ...
    'plotting', 1, 'hp_filt', [.4, .2, 3, 30], 'lp_filt', [17 25 3 50], ...
    'disp_ival', [-200 800],'varReject_freq_band', [5 25],  'ref_ival', []);

stdERPplots(epo, epo_r)

epo = proc_selectChannels(epo, scalpChannels);

epo_avg = proc_average(epo);


% data of non-target
dmy = proc_selectClasses(epo_avg, 'Non-Target');
% ERP-DATA in dmy.x
