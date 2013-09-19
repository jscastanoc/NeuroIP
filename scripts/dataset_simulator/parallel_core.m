function  parallel_core(model, dir_base, options)

nip_init();
% try
rng('default')
rng('shuffle')

t = model.t;
L = model.L;
fs = model.fs;
cortex = model.cortex;
clear model;

if isfield(cortex, 'vc') && isfield(cortex,'tri')
    cortex.vertices = cortex.vc;
    cortex.faces = cortex.tri;
end

Ntrials = options.Ntrials;
phase_shift = options.phase_shift;
snr_bio = options.snr_bio;
n_exp = options.n_exp;
Nspurious =options.Nspurious;
snr_meas = options.snr_meas;

for k = 1:length(snr_bio)
    for j = 1:n_exp
        cur_jobs = [];
        copy_res = {};
        
        % Phase shift for the sources (in seconds)
        ps= phase_shift + 0.05*randn(1,length(phase_shift));
        Nact = length(phase_shift);
        
        % Central frequency for the wavelet
        fc_wl = 9*ones(1,length(ps)) + 2.5*randn(1,length(ps));
        
        source_act = [];
        for m = 1:length(ps)
            source_act = [source_act ;aux_wavelet(fc_wl(m),t,ps(m))];
        end
        for i = Ntrials
            [y, Jclean, actidx] = nip_simtrials(L, cortex.vc, ...
                source_act, t, Nspurious , i, snr_meas, snr_bio(k),{'sample_all',false});
            dir = strcat(dir_base,num2str(length(ps)));
            file_name = strcat(dir,'/Exp',num2str(j),'Ntrials',...
                num2str(i),'BioNoise',num2str(snr_bio(k)),'.mat');
            save(file_name,'y','Jclean','actidx');
        end
    end
end
%     success = true;
% catch
%     success = false;
% end
end