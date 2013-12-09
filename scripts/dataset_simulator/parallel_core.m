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


[Nc, Nd] = size(L);
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

for kk = 1:n_exp
    % Phase shift for the sources (in seconds)
    ps= phase_shift + 0.05*randn(1,length(phase_shift));
    Nact = length(phase_shift);
    
    % Central frequency for the wavelet
    fc_wl = 9*ones(1,length(ps)) + 2.5*randn(1,length(ps));
    
    source_act = [];
    for m = 1:length(ps)
        source_act = [source_act ;aux_wavelet(fc_wl(m),t,ps(m))];
    end
    
    dir = randn(Nact,3);
    Nact = size(source_act,1);
    [Jclean, actidx] = nip_simulate_activity(cortex.vertices, Nact, source_act, dir, t,options);
    
    fuzzy = nip_fuzzy_sources(cortex,1.5);
    index = (1:3:Nd);
    for i = 0:2
        J(index+i,:) = fuzzy*Jclean(index+i,:); % J simulado FINAL
    end
    Jnorm = sqrt(sum(J.^2,2));
    idxaux = find(Jnorm<0.05*max(Jnorm));
    J(idxaux,:) = 0.0;
    CleanBrain{kk}.Jclean = sparse(J);
    CleanBrain{kk}.actidx = actidx;
end

for k = 1:length(snr_bio)
    for j = 1:n_exp
        cur_jobs = [];
        copy_res = {};
        for i = Ntrials
            [y, Jclean, actidx] = nip_simtrials(L, cortex.vc, CleanBrain{j}, t, Nspurious , i, snr_meas, snr_bio(k),{'sample_all',false});
            dir = strcat(dir_base,num2str(length(ps)));
            file_name = strcat(dir,'/Exp',num2str(j),'Ntrials',...
                num2str(i),'BioNoise',num2str(snr_bio(k)),'.mat');
            gof = norm(y-L*Jclean, 'fro')/norm(y, 'fro');
            save(file_name,'y','Jclean','actidx','fs','gof');
        end
    end
end
%     success = true;
% catch
%     success = false;
% end
end
