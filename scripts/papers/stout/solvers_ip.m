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
        fuzzy = nip_fuzzy_sources(model.cortex,i,struct('database','montreal','save',true));
        basisn = fuzzy(:,randi([1,model.Nd/3],nbasis,1));
        basisn = basisn/norm(basisn(:),1);
        basis = [basis basisn];
        group = [group n*ones(1,nbasis)];
        n = n+1;
    end
end

if sum(ismember(method,{'STOUT','TF-MxNE'}))
%     ltfatstart;
end

switch method
    case 'LOR'
        Q = speye(model.Nd);
        [J_rec,~] = nip_loreta(model.y,model.L,Q);        
    case 'S-FLEX'
       	[J_rec,~]  = nip_sflex(model.y,model.L,basis,5e-5);
    case 'TF-MxNE'
        [J_rec,~] = nip_tfmxne_port(model.y,model.L,[]);        
    case 'STOUT'
        [J_rec,~] = nip_stout(model.y,model.L,basis,[]);
    otherwise
        error(strcat('Ups! ',method,' is not available'))
end
Nd= size(model.cortex.vc,1);
distmat = graphrbf(model.cortex);
er = [];
time = [];
return

y = model.y;
L = model.L;
cortex = model.cortex;
clear model;
sig1 = sqrt(sum(J_rec.^2,2));
sig1 = sig1/norm(sig1);
sig2 = sqrt(sum(Jclean.^2,2));
sig2 = sig2/norm(sig2);
er(1) = nip_emd(sig1,sig2,distmat);
Jrecbckp = J_rec;
J_rec = J_rec/max(abs(J_rec(:)));
Jclean = Jclean/max(abs(Jclean(:)));
for iact = 1:actsources
    idxs = ((actidx(iact)-1)*3:(actidx(iact)-1)*3+2)+1;
    simact = repmat(Jclean(idxs,:),Nd,1);
    corr = sum(J_rec.*simact,2);
    corr = mean(reshape(corr,3,[]),1);
    [cormax(iact), idx_act(iact)] = max(corr);
    distact(iact) = distmat(idx_act(iact),actidx(iact));
end
er(2) = mean(cormax);
er(3) = mean(distact);
er(4) = mean((1./cormax).*distact);
er(5) = nip_error_tai(y,L,Jrecbckp);
[er(6),~,~] = nip_error_sai(cortex, Jclean,J_rec,5);

J_rec = sparse(J_rec);
time = toc;
end