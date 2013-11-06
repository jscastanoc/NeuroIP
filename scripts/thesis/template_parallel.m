% This is a template to run n tasks in parallel in the local computer (user
% chooses n)
% This template is to be used in my thesis.

clear; clc; close all;
<<<<<<< HEAD
nip_init
Nexp = 47:50 ;
% Nexp = 1;
Nact = [1,3,5];
% Nact = 3;
% Ntrials = [5 20 50 100 250];
Ntrials = [5 20 50 100 250];

% Ntrials = 250;
methods = {'KAL'};
snr_bio = -5;
dir_data = '/mnt/data/Master_Results/Datasets/simulated/montreal_sampleall_false/';
dir_results ='/mnt/data/Master_Results/TV-PRIORS/bin/montreal_sampleall_false/';
dir_error = '/mnt/data/Master_Results/TV-PRIORS/error/montreal_sampleall_false/';
Nparallel = 1;
=======

Nexp = 1:5 ;
Nact = [1,3];
Ntrials = [5];
methods = {'LOR','S-FLEX'};


icat = [];
jcat = [];
kcat = [];
mcat = [];

Nparallel = 7;
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
n = 0;
sched = findResource('scheduler','configuration','local');
job = createJob(sched);
compute = false;
for i = Nexp
<<<<<<< HEAD
    i
=======
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
    %     icat = cat(1,icat,i);
    for j = Nact
        %         jcat = cat(1,kcat,j);
        for k = Ntrials
            %             kcat = cat(1,kcat,k);
            for m = 1:numel(methods)
                %                 mcat = cat(1,mcat,m);
                n = n+1;
<<<<<<< HEAD
                args = {i,j,k,methods(m),dir_data,dir_results,dir_error,snr_bio};
                createTask(job, @thesis_core, 0,args);
=======
                args = {i,j,k,methods(m)};
                createTask(job, @dummy_core, 0,args);
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
                
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
