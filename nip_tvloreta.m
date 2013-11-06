function [J_est, extras] = nip_tvloreta(y,t,L,basis,cortex)
<<<<<<< HEAD
% [J_est, extras] = nip_tvloreta(y,t,L,basis,cortex)
=======
% J_est = nip_tvloreta(y,L,basis)
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e

[Nc,Nd] = size(L);
Nt = size(y,2);


<<<<<<< HEAD
[y_proj, y_rec, Ur, Er] = nip_tempcomp(y, t, [], 0.9);
Np = size(Ur,2);


aff = nip_fuzzy_sources(cortex,4,struct('dataset','montreal','save',true));

finalD =[];
for i = 1:Np
    %     D(:,i) = nip_lcmv(y_proj(:,i),L);
    D(:,i) = nip_sflex(y_proj(:,i),L,basis, 60);
    
    if isempty(find(D(:,i)))
        D(:,i) = [];
    else
        % Energy Threshold
        De = nip_energy(D(:,i));
        idx = find(De < 0.01*max(De));
        De(idx) = 0;
        for j = 1:3
            idx3 = (idx-1)*3+j;
            D(idx3,i) = 0;
        end
        D(:,i) = D(:,i)./norm(D(:,i));
        %
        idx = find(De);
        affcr = aff(idx,idx);
        S = svd(affcr);
        cum_var = cumsum(diag(S))/sum(diag(S));
        Nclusters = numel(find(cum_var <= 0.5));
        if Nclusters == 0;
            Nclusters = 1;
        end
        label = knkmeans(affcr,Nclusters);
        
        
        labelfull = zeros(Nd/3,1);
        labelfull(idx) = label;
        
        finalDi = zeros(size(D(:,i),1),Nclusters+1);
        
        n = 1;
        for j = unique(labelfull)'
            idx1 = find(labelfull == j);
            for k = 1:3
                idx3 = (idx1-1)*3+k;
                finalDi(idx3,n) = D(idx3,i);
            end
            n = n+1;
        end
        finalD = [finalD,finalDi(:,2:end)];
    end
=======
[y_proj, y_rec, Ur, Er] = nip_tempcomp(y, t, [0.1 40], 2);
Np = size(Ur,2);


aff = nip_fuzzy_sources(cortex,5,struct('dataset','montreal','save',true));

finalD =[];
for i = 1:Np
%     D(:,i) = nip_lcmv(y_proj(:,i),L);
    D(:,i) = nip_sflex(y_proj(:,i),L,basis, 1e-4);
    
    % Energy Threshold
    De = nip_energy(D(:,i));
    idx = find(De < 0.1*max(De));
    De(idx) = 0;
    for j = 1:3
        idx3 = (idx-1)*3+j;
        D(idx3,i) = 0;
    end
    D(:,i) = D(:,i)./norm(D(:,i));
    
    idx = find(De);
    affcr = aff(idx,idx);
    S = svd(affcr);
    
    Nclusters = numel(find(S >= 0.6*max(S(:))));
    label = knkmeans(affcr,Nclusters);
    
    
    labelfull = zeros(Nd/3,1);
    labelfull(idx) = label;
    
    finalDi = zeros(size(D(:,i),1),Nclusters+1);
    
    n = 1;
    for j = unique(labelfull)'
        idx1 = find(labelfull == j);
        for k = 1:3
            idx3 = (idx1-1)*3+k;
            finalDi(idx3,n) = D(idx3,i);
        end
        n = n+1;
    end
    finalD = [finalD,finalDi(:,2:end)];
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
end

for i=1:size(finalD,2)
    finalD(:,i) = finalD(:,i)./norm(finalD(:,i));
end
<<<<<<< HEAD
finalD = D;
=======
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e

idx = 1:3:Nd;
n = 1;
clear D
D = zeros(Nd,size(finalD,2)*3);
<<<<<<< HEAD
for i = 1:size(finalD,2);
    fDe(:,i) = nip_energy(finalD(:,i));
    for j = 0:2
        D(idx+j,n) = fDe(:,i);
        n = n+1;
=======
for i = 1:size(finalD,2);  
    fDe(:,i) = nip_energy(finalD(:,i));
    for j = 0:2
            D(idx+j,n) = fDe(:,i);
            n = n+1;
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
    end
end

% Matrix with the perfectly located dictionaries mapped in to the lead
% field
Laux = L*D;
h = nip_kalman_hyper(y,Laux);

for i = 1:size(Ur,2)
    hp = h*Ur(:,i);
    Q = [];
    for j = 1:size(D,2)
<<<<<<< HEAD
        Q = [Q, abs(hp(j))*D(:,j)];
    end
    Q = sum(Q,2);
    %     Q = ones(size(Q));
=======
        Q = [Q, hp(j)*D(:,j)];
    end
    Q = sum(Q,2);
    
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
    [J_est(:,i),~] = nip_loreta(y_proj(:,i),L,sparse(diag(Q)));
end

J_est = J_est*Ur';
extras.D = D;
extras.h = h;