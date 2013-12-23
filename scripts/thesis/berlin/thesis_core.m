function [J_rec, time, er, extra] = thesis_core(Nexp,Nact,Ntrials,methods,dir_data,dir_results,dir_error,snr_bio,model, Jclean, actidx,depth)

nip_init();
% load_data;




switch depth
    case 'none'
        L = model.L;
    case 'Lnorm'
        gamma = 0.6;
        [L, extras] = nip_depthcomp(model.L,struct('type',depth,'gamma',gamma));
        Winv = extras.Winv;
    case 'sLORETA'
        [L, extras] = nip_depthcomp(model.L,struct('type',depth));
        Winv = extras.Winv;
end
clear extras;

for m = 1:numel(methods)
    tic;
    switch methods{m}
        case 'LOR'
            %                         Q = nip_translf(repmat(nip_translf(repmat(basis,[1 1 3]))',[1 1 3]));
            Q = speye(size(L,2));
            [J_rec,extra]  = nip_loreta(model.y,L,Q);
        case 'KAL'
            %                         neigh = nip_translf(repmat(nip_translf(repmat(basis,[1 1 3]))',[1 1 3]));
            par = [0.5,0.01,0.002,1e-4,1e-2];
            neigh = speye(model.Nd);
            [J_rec,extra]  = nip_kalmanwh(model.y,L,neigh,par);
        case 'IRA3'
            neigh = speye(model.Nd);
            Q = speye(model.Nd);
            pariter = [0.0013 0.0084];
            [J_rec,extra] = nip_iterreg(model.y,model.L,Q,neigh,3,eye(model.Nc),pariter);
        case 'IRA5'
            neigh = speye(model.Nd);
            Q = speye(model.Nd);
            pariter = [0.0013 0.0084];
            [J_rec,extra] = nip_iterreg(model.y,model.L,Q,neigh,5,eye(model.Nc),pariter);
        case 'LOR_PROJ'
            %                         Q = nip_translf(repmat(nip_translf(repmat(basis,[1 1 3]))',[1 1 3]));
            Q = speye(model.Nd);
            [y_proj,~,Ur,~]= nip_tempcomp(model.y,model.t,[0 60],0.9);
            [J_rec,extra] = nip_loreta(y_proj,model.L,Q);
            J_rec = Jrec*Ur';
        case 'LORTV'
            [J_rec,extra]=nip_tvloreta(model.y,model.t,model.L,basis,model.cortex);
        case 'S-FLEX'
            % Spatial dictionary
            sigma = 1;
            B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));
            B = nip_blobnorm(B,'norm',2);
            % Options for the inversion
            reg_par = 100;
            [J_rec, extras] = nip_sflex(model.y, L, B, 'regpar', reg_par,'optimgof',true,'gof',0.8);
        case 'TF-MxNE'
            spatial_reg = 100;
            temp_reg =  5;
            [J_rec, extras] = nip_tfmxne_port(model.y, L, 'optimgof',true,...
                'sreg',spatial_reg,'treg',temp_reg,'gof',0.9);
        case 'STOUT'
            % Spatial dictionary
            sigma = 1;
            B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));
            %         n = 10;
            %         idx = randsample(4001,1000);
            B = nip_blobnorm(B,'norm',2);
            
            
            % Options for the inversion
            spatial_reg = 100;
            temp_reg =  5;
            [J_rec, extras] = nip_stout(model.y, L, B,'optimgof',false,...
                'sreg',spatial_reg,'treg',temp_reg,'gof',0.9);
        otherwise
            error('%s not available as method',methods{j})
    end
    idx = find(sqrt(sum(J_rec.^2)) <= 0.01*sqrt(sum(J_rec.^2)));
    J_rec(idx,:) = zeros(length(idx),model.Nt);
    
    
    if ~strcmp(depth,'none')
        J_est = nip_translf(J_rec');
        for i = 1:3
            J_est(:,:,i) = (Winv(:,:,i)*J_est(:,:,i)')';
        end
        J_rec = nip_translf(J_est)';
    end
    
    J_rec = full(J_rec);
    Jclean = full(Jclean);
    %%% Compute Errors %%%
    er = nip_all_errors(model.y(:,round(end/9):end),model.L,...
        J_rec(:,round(end/9):end),Jclean(:,round(end/9):end),model.cortex,actidx);
    %%% -------------- %%%
    
    
    J_rec = sparse(J_rec);
    time = toc;
    %                 dir = strcat(dir_results,num2str(j));
    %                 file_name = strcat(dir,'/',methods{m},'Exp',num2str(i),'Ntrials',...
    %                     num2str(k),'BioNoise',num2str(snr_bio),'.mat');
    %
    %                 save(file_name,'J_rec','extra','time');
    %                 dir = strcat(dir_error,num2str(j));
    %                 file_name = strcat(dir,'/',methods{m},'Exp',num2str(i),'Ntrials',...
    %                     num2str(k),'BioNoise',num2str(snr_bio),'.mat');
    %                 save(file_name,'er');
end
end
