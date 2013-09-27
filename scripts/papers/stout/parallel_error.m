function  parallel_error(method, dir_sim, dir_results, dir_error, options)

nip_init();
% try
rng('default')
rng('shuffle')

cortex = options.cortex;
clear model;

if isfield(cortex, 'vc') && isfield(cortex,'tri')
    cortex.vertices = cortex.vc;
    cortex.faces = cortex.tri;
end

Ntrials = options.Ntrials;
snr_bio = options.snr_bio;
n_exp = options.n_exp;
Nspurious =options.Nspurious;
snr_meas = options.snr_meas;
L = options.L;
clear options 

Nd = size(cortex.vertices,1);
distmat = graphrbf(cortex);
act_sources = [1,3,5];
for j = 1:n_exp
    cur_jobs = [];
    copy_res = {};
    for l = 1:length(act_sources)
        for i = Ntrials
            dir = strcat(dir_results,num2str(act_sources(l)));
            file_name = strcat(dir,'/',method,'Exp',num2str(j),'Ntrials',...
                num2str(i),'BioNoise',num2str(snr_bio),'.mat');
            load(file_name)
            dir = strcat(dir_sim, num2str(act_sources(l)));
            file_name = strcat(dir,'/Exp',num2str(j),'Ntrials',...
                num2str(i),'BioNoise',num2str(snr_bio),'.mat');
            load(file_name);
            
            sig1 = sqrt(sum(J_rec.^2,2));
            sig1 = sig1/norm(sig1);
            sig2 = sqrt(sum(Jclean.^2,2));
            sig2 = sig2/norm(sig2);
            er(1) = nip_emd(sig1,sig2,distmat);
            Jrecbckp = J_rec;
            J_rec = J_rec/max(abs(J_rec(:)));
            Jclean = Jclean/max(abs(Jclean(:)));
            for iact = 1:length(act_sources(l))
                idxs = ((actidx(iact)-1)*3:(actidx(iact)-1)*3+2)+1;
                simact = repmat(Jclean(idxs,:),Nd,1);
                corr = sum(J_rec.*simact,2);
                corr = mean(reshape(corr,3,[]),1);
                [cormax(iact), idx_act(iact)] = max(corr);
                distact(iact) = distmat(idx_act(iact),actidx(iact));
            end
            er(2) = mean(cormax);
            er(3) = mean(distact);
            er(4) = mean((1/cormax).*distact);
            er(5) = nip_error_tai(y,L,Jrecbckp);
            
            dir = strcat(dir_error,num2str(act_sources(l)));
            file_name = strcat(dir,'/',method,'Exp',num2str(j),'Ntrials',...
                num2str(i),'BioNoise',num2str(snr_bio),'.mat');
            save(file_name,'er');
        end
    end
    
end
end