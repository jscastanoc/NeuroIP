<<<<<<< HEAD
function run_parallel_exp(Nexp,Nact,Ntrials,methods,Nparallel, dir_data,dir_results, dir_error, snr_bio)
=======
function run_parallel_exp(Nexp,Nact,Ntrials,methods,Nparallel, dir_data,dir_results, snr_bio)
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
% Function to run in parallel the experiments using (already existent)
% simulated data corresponding to:
% Input:
%       Nexp -> Run experiments number Nexp(1)...Nexp(N)
%       Nact -> Run experiments for Nact(1)...Nact(M) active sources
%       Ntrials -> Run experiments for Ntrials(1)...Ntrials(P)
%       methods -> cell containing the names of the methods to use as solvers
%       snr_bio -> Biological noise added to the simulation
%       dir_data -> directory where the simulations are
%       res_data -> directory where the results are going to be saved



n = 0;
sched = findResource('scheduler', 'configuration', 'local');
job = createJob(sched);
for i = Nexp
    for j = Nact
        for k = Ntrials
            for m = 1:numel(methods)
                n = n+1;
<<<<<<< HEAD
                args = {i,j,k,methods(m),dir_data,dir_results,dir_error,snr_bio};
=======
                args = {i,j,k,methods(m),dir_data,dir_results,snr_bio};
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
                createTask(job, @core, 0,args);
                
                if ((n >= Nparallel) || ((i == Nexp(end)) && ...
                        (j == Nact(end)) && (k == Ntrials(end)) && m == numel(methods)))
                    fprintf('Computing %u tasks in parallel\n',n);
                    fprintf('Nexp: %u',i)
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
end

<<<<<<< HEAD
function core(Nexp,Nact,Ntrials,methods,dir_data,dir_results,dir_error,snr_bio)
=======
function core(Nexp,Nact,Ntrials,methods,dir_data,dir_results,snr_bio)
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
nip_init();
load_data;

% Depth compensation
<<<<<<< HEAD
[model.L, extras] = nip_depthcomp(model.L, struct('type','sLORETA') );
% [model.L] = nip_depthcomp(model.L, struct('type','Lnorm','gamma',0.8) );
Winv =extras.Winv;
clear extras;
=======
% [model.L, extras] = nip_depthcomp(model.L, struct('type','sLORETA') );
% [model.L] = nip_depthcomp(model.L, struct('type','Lnorm','gamma',0.8) );
% Winv =extras.Winv;
% clear extras;
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e

if sum(ismember(methods,{'LOR','KAL','IRA3','IRA5','LORPROJ'}))
    Laplacian = laplace_loreta(model.cortex.vc);
    Lapaux = Laplacian;
    Laplacian = repmat(full(Laplacian),[1,1,3]);
    Laplacian = nip_translf(Laplacian);
    Laplacian = repmat(Laplacian',[1,1,3]);
    Laplacian = sparse(nip_translf(Laplacian));
    Q = speye(model.Nd);
end

if sum(ismember(methods,{'STOUT','SFLEX'}))
    basis = nip_fuzzy_sources(model.cortex,2,struct('dataset','montreal','save',true));
end

if sum(ismember(methods,{'TF-MxNE','STOUT'}))
    ltfatstart;
end


for i = Nexp'
    for j = Nact'
        for k = Ntrials'
            dir = strcat(dir_data,num2str(j));
            file_name = strcat(dir,'/Exp',num2str(i),'Ntrials',...
                num2str(k),'BioNoise',num2str(snr_bio),'.mat');
            load(file_name)
            model.fs = fs;
            model.y = y;
            model.Nt = size(y,2);
            model.t = 0:1/fs:model.Nt/fs;
            for m = 1:numel(methods)
                tic;
                switch methods{m}
                    case 'LOR'
                        [J_est,extra]  = nip_loreta(model.y,model.L,Q);
                    case 'KAL'
                        %                         par =
                        [J_est,extra]  = nip_kalmanwh(model.y,model.L,Laplacian,par);
                    case 'IRA3'
                        %                         par
                        [J_est,extra] = nip_iterreg(model.y,model.L,Q,Laplacian,3,eye(model.Nc),par);
                    case 'IRA5'
                        %                         par
                        [J_est,extra] = nip_iterreg(model.y,model.L,Q,Laplacian,3,eye(model.Nc),par);
                    case 'LOR_PROJ'
                        [y_proj,~,Ur,~]= nip_tempcomp(model.y,model.t,[0 60],0.9);
                        [J_est,extra] = nip_loreta(y_proj,model.L,Q);
                        J_est = Jrec*Ur';
                    case 'S-FLEX'
                        [J_est, extra] = nip_sflex(model.y, model.L, basis, reg_par);
                    case 'TF-MxNE'
<<<<<<< HEAD
                        [J_est,extra] = nip_tfmxne_port(model.y, model.L, options,[]);
=======
                        [J_est,extra] = nip_tfmxne_port(model.y, model.L, options);
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
                    case 'STOUT'
                        [J_est,extra] = nip_stout(model.y, model.L, basis,[]);
                    otherwise
                        err('%s not available as method',methods{j})
                end
                time = toc;
                dir = strcat(dir_results,num2str(j));
                file_name = strcat(dir,'/',methods{m},'Exp',num2str(i),'Ntrials',...
                    num2str(k),'BioNoise',num2str(snr_bio),'.mat');
                
<<<<<<< HEAD
                % ----- Uncomment if depth compensation with sLORETA ----%
                J_est = nip_translf(J_est');
                J_est = permute(J_est,[2 1 3]);
                for i = 1:3
                    J_est(:,:,i) = Winv(:,:,i)*J_est(:,:,i);
                end
                J_est = permute(J_est,[2 1 3]);
                J_est = nip_translf(J_est)';
=======
                % -------Comment if depth compensation with sLORETA------%
%                 J_est = nip_translf(J_est');
%                 J_est = permute(J_est,[2 1 3]);
%                 for i = 1:3
%                     J_est(:,:,i) = Winv(:,:,i)*J_est(:,:,i);
%                 end
%                 J_est = permute(J_est,[2 1 3]);
%                 J_est = nip_translf(J_est)';
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
                %--------------------------------------------------------%
                
                idx = find(sqrt(sum(J_est.^2)) <= 0.01*sqrt(sum(J_est.^2)));
                J_est(idx,:) = zeros(length(idx),model.Nt);
                J_est = sparse(J_est);
                
                save(file_name,'J_est','extra','time');
<<<<<<< HEAD
                
                
                dir = strcat(dir_error,num2str(j));
                file_name = strcat(dir,'/',methods{m},'Exp',num2str(i),'Ntrials',...
                    num2str(k),'BioNoise',num2str(snr_bio),'.mat');
                er = nip_all_errors(model.y,model.L,J_est,Jclean,model.cortex,actidx);
                save(file_name,'er');
=======
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
            end
        end
    end
end

end