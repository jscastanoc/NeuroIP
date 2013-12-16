
function run_parallel_exp(Nexp,Nact,Ntrials,methods,Nparallel, dir_data,dir_results, dir_error, snr_bio,depth)
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
                args = {i,j,k,methods(m),dir_data,dir_results,dir_error,snr_bio,depth};

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

function core(Nexp,Nact,Ntrials,methods,dir_data,dir_results,dir_error,snr_bio,depth)

nip_init();
load_data;

switch depth
    case 'none'
        L = model.L;
    case 'Lnorm'        
        gamma = 0.8;
        [L, extras] = nip_depthcomp(model.L,struct('type',depth,'gamma',gamma));
        Winv = extras.Winv;
    case 'sLORETA'
        [L, extras] = nip_depthcomp(model.L,struct('type',depth));
        Winv = extras.Winv;        
end
clear extras;


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
						Q = eye(model.Nd); %Matriz de covarianza apriori
						[J_est, extras] = nip_loreta(model.y, L, Q);
					case 'TF-MxNE'
						spatial_reg = 100;
						temp_reg =  1;
						options.iter = 50;
						options.tol = 2e-2;
						[J_est, extras] = nip_tfmxne_port(model.y, L, 'optimgof',true,...
								'sreg',spatial_reg,'treg',temp_reg);
					case 'STOUT' 
						% Spatial dictionary
						sigma = 1;
						B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));
				%         n = 10;
				%         idx = randsample(4001,1000);        
						B = nip_blobnorm(B);
						
						
						
						% Options for the inversion
						spatial_reg = 250;
						temp_reg =  1;
						options.iter = 50;
						options.tol = 2e-2;
						[J_est, extras] = nip_stout(model.y, L, B,'optimgof',true,...
								'sreg',spatial_reg,'treg',temp_reg);
					case 'S-FLEX'
						% Spatial dictionary
						sigma = 1;        
						B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));
						B = nip_blobnorm(B);
						% Options for the inversion
						reg_par = 100;
						[J_est, extras] = nip_sflex(model.y, L, B, 'regpar', reg_par,'optimgof',true);
				end
				
                time = toc;
                dir = strcat(dir_results,num2str(j));
                file_name = strcat(dir,'/',methods{m},'Exp',num2str(i),'Ntrials',...
                    num2str(k),'BioNoise',num2str(snr_bio),'.mat');
                

               
				if ~strcmp(depth,'none')
					J_est = nip_translf(J_est');
					for i = 1:3
						J_est(:,:,i) = (Winv(:,:,i)*J_est(:,:,i)')';
					end
					J_est = nip_translf(J_est)';
				end


                %--------------------------------------------------------%
                
                idx = find(sqrt(sum(J_est.^2)) <= 0.01*sqrt(sum(J_est.^2)));
                J_est(idx,:) = zeros(length(idx),model.Nt);
                J_est = sparse(J_est);
                
                save(file_name,'J_est','extra','time');
           
                
                dir = strcat(dir_error,num2str(j));
                file_name = strcat(dir,'/',methods{m},'Exp',num2str(i),'Ntrials',...
                    num2str(k),'BioNoise',num2str(snr_bio),'.mat');
                er = nip_all_errors(model.y,model.L,J_est,Jclean,model.cortex,actidx);
                save(file_name,'er');

            end
        end
    end
end

end
