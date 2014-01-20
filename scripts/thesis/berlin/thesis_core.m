function [J_rec, time, er, extras] = thesis_core(Nexp,Nact,Ntrials,method,dir_data,dir_results,dir_error,snr_bio,model, Jclean, actidx,depth, resnorm)

nip_init();
% load_data;




switch depth
    case 'none'
        L = model.L;
        Winv =[];
    case 'Lnorm'
        gamma = 0.6;
        [L, extras] = nip_depthcomp(model.L,struct('type',depth,'gamma',gamma));
        Winv = extras.Winv;
    case 'sLORETA'
        [L, extras] = nip_depthcomp(model.L,struct('type',depth));
        Winv = extras.Winv;
end
clear extras;

tic;
switch method
    case 'LOR'
        [J_rec, extras] = nip_loreta(model.y,L,'cov',speye(size(L,2)),'Winv',Winv);
    case 'S-FLEX'
        % Spatial dictionary
        sigma = 1;
        B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));
        B = nip_blobnorm(B,'norm',2);
        % Options for the inversion
        reg_par = 100;
        [J_rec, extras] = nip_sflex(model.y, L, B,'optimres',true,'regpar', reg_par,'resnorm',resnorm,'Winv',Winv);
    case 'TF-MxNE'
        spatial_reg = 160;
        temp_reg =  10;
        [J_rec, extras] = nip_tfmxne_python(model.y, L, 'optimres',true,'sreg',spatial_reg,'treg',temp_reg,'resnorm',resnorm,'Winv',Winv);
    case 'STOUT'
        % Spatial dictionary
        sigma = 1;
        B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));
        B = nip_blobnorm(B,'norm',2);
        
        % Options for the inversion
        spatial_reg = 160;
        temp_reg =  10;
        [J_rec, extras] = nip_stout_python(model.y, L, B,'optimres',true,'sreg',spatial_reg,'treg',temp_reg,'resnorm',resnorm,'Winv',Winv);
    otherwise
        error('%s not available as method',methods{m})
end

J_rec = full(J_rec);
Jclean = full(Jclean);

J_rec = J_rec*(norm(model.y,'fro')/norm(model.L*J_rec,'fro'));
extras.resnorm = norm(model.y-model.L*J_rec, 'fro')/norm(model.y, 'fro');
extras.reser = resnorm-extras.resnorm;
%%% Compute Errors %%%
er = nip_all_errors(model.y(:,round(end/9):end),model.L,...
    J_rec(:,round(end/9):end),Jclean(:,round(end/9):end),model.cortex,actidx);

J_rec = sparse(J_rec);
time = toc;
end
