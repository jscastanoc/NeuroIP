function [J_rec time err] = core_matgrid_test(model,method, Jclean, actsources, actidx)

% Width of the basis functions for S-FLEX and stout
iter_basis = [1.5];

if sum(ismember(method,{'S+T','S-FLEX'}))
    nbasis = size(model.L,2)/3;
    basis = [];
    n = 1;
    group = [];
    for i = iter_basis
        fuzzy = nip_fuzzy_sources(model.cortex,i);
        basisn = fuzzy(:,randi([1,model.Nd/3],nbasis,1));
        basisn = basisn/norm(basisn(:),1);
        basis = [basis basisn];
        group = [group n*ones(1,nbasis)];
        n = n+1;
    end
end

switch method
    case 'LOR'
        %         [Laplacian] = nip_neighbor_mat(model.cortex);
        %         Q = inv(Laplacian*Laplacian');
        Q = speye(model.Nd);
        [J_rec,~] = nip_loreta(model.y,model.L,Q);
        
    case 'S-FLEX'
        index = (1:3:model.Nd);
        for i = 0:2
            L(:,index+i) = model.L(:,index+i)*basis;
        end
        [xx,~] = sflex_cortical_dal(model.y,L,[],struct('eps',0.25,'B',eye(model.Nd)));
        J_rec = permute(xx,[1 3,2]);
        index = (1:3:model.Nd);
        for i = 0:2
            J_rec(index+i,:) = basis*J_rec(index+i,:);
        end
        
    case 'TF-MxNE'
        options.iter = 50;
        options.spatial_reg = 0.9;
        options.temp_reg = 0.15;
        options.tol = 1e-2;
        options.a = 10;
        options.M = 200;
        [J_rec,~] = nip_tfmxne_port(model.y,model.L,options);
        
    case 'S+T'
        options.iter = 50;
        options.spatial_reg = 0.9;
        options.temp_reg = 0.15;
        options.tol = 1e-2;
        options.a = 10;
        options.M = 200;
        [J_rec,~] = nip_sflex_tfmxne(model.y,model.L,basis,options);
    otherwise
        error(strcat('Nah! ',method,' is not available'))
end
Nd = size(model.cortex.vc,1);
distmat = graphrbf(model.cortex);

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
er(4) = mean((1/cormax).*distact);
er(5) = nip_error_tai(y,L,Jrecbckp);

J_rec = sparse(J_rec);
time = toc;
end