function [J_rec time er] = solvers_ip(model,method, Jclean, actsources, actidx)

nip_init();
% Width of the basis functions for S-FLEX and stout
iter_basis = [1];

if sum(ismember(method,{'STOUT','S-FLEX'}))
    nbasis = size(model.L,2)/3;
    basis = [];
    n = 1;
    group = [];
    for i = iter_basis
        fuzzy = nip_fuzzy_sources(model.cortex,i,struct('dataset','montreal','save',true));
        basisn = fuzzy(:,randi([1,model.Nd/3],nbasis,1));
        basisn = basisn/norm(basisn(:),1);
        basis = [basis basisn];
        group = [group n*ones(1,nbasis)];
        n = n+1;
    end
end

if sum(ismember(method,{'STOUT','TF-MxNE'}))
    ltfatstart;
end

tic

switch method
    case 'LOR'
        Q = eye(model.Nd); %Matriz de covarianza apriori
        [J_est, extras] = nip_loreta(model.y, L, Q);
    case 'TF-MxNE'
        spatial_reg = 250;
        temp_reg =  1;
        options.iter = 50;
        options.tol = 2e-2;
        [J_est, extras] = nip_tfmxne_port(model.y, L, 'optimgof',true,...
            'sreg',spatial_reg,'treg',temp_reg,'gof',model.gof);
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
            'sreg',spatial_reg,'treg',temp_reg,'gof',model.gof);
    case 'S-FLEX'
        % Spatial dictionary
        sigma = 1;
        B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));
        B = nip_blobnorm(B);
        % Options for the inversion
        reg_par = 100;
        [J_est, extras] = nip_sflex(model.y, L, B, 'regpar', reg_par,'optimgof',true,'gof',model.gof);
end
J_rec = J_est;
clear J_est;
Nd= size(model.cortex.vc,1);
distmat = graphrbf(model.cortex);
er = [];

time = toc;

y = model.y;
L = model.L;
cortex = model.cortex;

er = nip_all_errors(y,L,J_rec,Jclean,cortex,actidx);

J_rec = sparse(J_rec);
time = toc;
end