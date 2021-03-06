% This is a template to run n tasks in parallel in the local computer (user
% chooses n)

clear; clc; close all;
nip_init()
Nexp = [1] ;
% Nexp = 1;
Nact = [1];
% Nact = 3;
% Ntrials = [5 20 50 100 250];
Ntrials = [250];

% Ntrials = 250;
methods = {'KAL','IRA3','IRA5','LOR_PROJ','LOR','DOBERMAN'};
methods = {'IRA3'}
snr_bio = -5;
dir_data = '/mnt/data/Master_Results/Datasets/simulated/montreal_sampleall_false/';
dir_results ='/mnt/data/Master_Results/IRA/bin/montreal_sampleall_false/';
dir_error = '/mnt/data/Master_Results/IRA/error/montreal_sampleall_false/';
Nparallel = 2;
n = 0;
sched = findResource('scheduler','configuration','local');
job = createJob(sched);
compute = false;
depth = 'Lnorm';
for i = Nexp
    for j = Nact
        for k = Ntrials
            for m = 1:numel(methods)
                n = n+1;
                dir = strcat(dir_data,num2str(j),'/');
                    file_name = strcat(dir,'Exp',num2str(i),'Ntrials',...
                        num2str(k),'BioNoise',num2str(snr_bio),'.mat');
                    load(file_name);
                    load_data;
                    model.fs = fs;
                    model.y = y;
                    model.Nt = size(y,2);
                    model.t = 0:1/fs:model.Nt/fs;
                    
                    resnorm = norm(model.y - model.L*Jclean, 'fro')/norm(model.y, 'fro')

                args = {i,j,k,methods{m},dir_data,dir_results,dir_error,snr_bio,...
                    model,Jclean,actidx,depth,resnorm};
                createTask(job, @thesis_core, 0,args);
                if ((n >= Nparallel) || ((i == Nexp(end)) && ...
                        (j == Nact(end)) && (k == Ntrials(end)) && m == numel(methods)))
                    fprintf('Computing %u tasks in parallel\n',n);
                    submit(job);
                    waitForState(job, 'finished');
                    errmsgs = get(job.Tasks, {'ErrorMessage'});
                    nonempty = ~cellfun(@isempty, errmsgs);
                    celldisp(errmsgs(nonempty));
                    destroy(job);
                    job = createJob(sched);
                    
                    % Cleanup
                    compute = false;
                    n = 0;
                end
            end
        end
    end
    
end
