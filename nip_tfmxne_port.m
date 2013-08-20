function [J_rec extras] = nip_tfmxne_port(y,L,options)
% function [J_rec extras] = nip_tfmxne_port(y,L,options)
% Solve the inverse problem using the TF-MxNE approach as introduced by
% Gramfort et al 2012.
%  Input:
%         y -> NcxNt. Matrix containing the data,
%         L -> Ncx3Nd. Lead Field matrix
%         options -> struct.
%                 options.spatial_reg -> scalar. Percentage of spatial
%                       regularization (between 0 and 1). 
%                 options.temp_reg -> scalar. Percentage of temporal
%                       regularization (between 0 and 1).
%                 options.a -> scalar. Time shift for the time frequency
%                       transform.
%                 options.m -> scalar. Frequency bins for the time frequency
%                       transform.
%                 options.tol -> Scalar. Default 1e-3
%   Output:
%         J_rec -> 3NdxNt. Reconstructed activity (solution)
%         extras. -> Currently empty
%
% Comments: 
% Almost a literal translation of the python function tf_mixed_norm_solver found in:
% https://github.com/mne-tools/mne-python/blob/0ac3ac1a1634673da013109f38daf4f162cae117/mne/inverse_sparse/mxne_optim.py
% Juan S. Castano C.
% jscastanoc@gmail.com
% 14 Aug 2013
 
[Nc Nd] = size(L);
Nt = size(y,2);

% Initialization of the TF-MxNE algorithm
if ~isfield(options,'a')||~isfield(options,'m')
    options.a = 10;
    options.m = 200;
end

if ~isfield(options,'spatial_reg');options.spatial_reg = 0.7;end
if ~isfield(options,'temp_reg');options.temp_reg = 0.05;end
if ~isfield(options,'iter');options.iter = 5;end
if ~isfield(options,'tol');options.tol = 1e-3;end

tol = options.tol;
spatial_reg = options.spatial_reg;
temp_reg = options.temp_reg;

a = options.a;
M = options.m;  

c = dgtreal(y','gauss',a,M);
T = size(c,2);
K = size(c,1);
Z = sparse(0,K*T);
Y = sparse(Nd,K*T);
J_rec = sparse(Nd,Nt);
J_recf = sparse(Nd,M);


tau = 1;
% lambda = 0.001; % Time regularization parameter
base_par = 0; % Spatial regularization parameter
tempGY = L'*y;
for i = 1:Nd
    basepar = max(norm(tempGY(i,:),2),base_par);
end
mu = basepar*spatial_reg; % Spatial regularization parameter
lambda = basepar*temp_reg; % Time regularization parameter
clear tempGY;

R = y;

active_set = logical(sparse(1,Nd));
Y_time_as = [];
Y_as = [];
lipschitz_k = lipschitz_contant(y, L, 1e-3, a, M);
mu_lc = mu/lipschitz_k;
lambda_lc = lambda/lipschitz_k;
stop = inf;
fprintf('Running TF-MxNE algorithm... \n');
rev_line = '';
eta = 0;
for i = 1:options.iter
    tic;
    msg = sprintf('Iteration # %d, Stop: %d, Elapsed time: %f',i,stop,eta);
    fprintf([rev_line, msg]);
    rev_line = repmat(sprintf('\b'),1,length(msg));
    
    Z_0 = Z;   active_set0 = active_set;
    
    if (sum(active_set) < size(R,1)) && ~isempty(Y_time_as)
        GTR = L'*R/lipschitz_k;
        A = GTR;
        A(find(Y_as),:) = A(find(Y_as),:) + Y_time_as(find(Y_as),1:Nt);
        [~, active_set_l21] = prox_l21(A,mu_lc,3);
        idx_actsetl21 = find(active_set_l21);
        
        aux = dgtreal(GTR(idx_actsetl21,:)','gauss',a,M);
        aux = permute(aux,[3 1 2]);
        aux = reshape(aux,sum(active_set_l21),[]);
        
        B = Y(idx_actsetl21,:) + aux;
        [Z, active_set_l1] = prox_l1(B,lambda_lc,3);
        active_set_l21(idx_actsetl21) = active_set_l1;
        active_set_l1 = active_set_l21;
    else
        temp = dgtreal(R','gauss',a,M);
        temp = permute(temp,[3 1 2]);
        temp = reshape(temp,Nc,[]);
        Y = Y + L'*temp/lipschitz_k;
        [Z, active_set_l1] = prox_l1(Y,lambda_lc,3);
    end
    [Z, active_set_l21] = prox_l21(Z,mu_lc,3);
    active_set = active_set_l1;
    active_set(find(active_set_l1)) = active_set_l21;
    
    
    
    if i < options.iter
        % line 7 of algorithm 1 (see Ref paper)
        tau_0 = tau;
        
        % line 8 of algorithm 1 (see Ref paper)
        tau = 0.5+sqrt(1+4*tau^2)/2;
        
        Y = sparse(Nd,K*T);
        
        dt = (tau_0-1)/tau;
        
        Y(find(active_set),:) = (1 + dt)*Z;
        Y(find(active_set0),:)= Y(find(active_set0),:) - dt*Z_0;

        
        Y_as = active_set0|active_set;

        temp =  reshape(full(Y)',K,T,[]);
        temp = flipdim(temp,2);
        Y_old = Y_time_as;
        Y_time_as = flipud(idgtreal(temp,'gauss',a,M))';
        
        if isempty(Y_time_as)||isempty(Y_old)
            stop_old = inf;
            stop = inf;
        else
            %        stop = (max(max(abs(Z)))-max(max(abs(Z_0))))/max(max(abs(Z_0)));
            %         stop = (norm(abs(Z))-norm(abs(Z_0)))/norm(abs(Z_0));
            stop_old = stop;
            stop = norm(Y_time_as - Y_old)/norm(Y_old);
        end        
        if stop < tol || stop_old < stop;
            Z = Z_0;
            active_set = active_set0;
            break
        end
        if stop < tol 
            break
        end
        
        R = y - L(:, find(Y_as))*Y_time_as(find(Y_as),1:Nt);
    end
    eta = toc;
end

fprintf(' \nDone!... \nTransforming solution to the time domain: \n%d non-zero time series \n'...
    , length(active_set))

temp =  reshape(full(Z)',K,T,[]);
temp = flipdim(temp,2);
J_recf(find(active_set),:) = flipud(idgtreal(temp,'gauss',a,M))';
J_rec = J_recf(:,1:Nt);

extras = [];
end

function [Y active_set] = prox_l21(Y,mu,n_orient)
    n_pos = size(Y,1)/n_orient;
    
    rows_norm = sqrt(sum(reshape(abs(Y).^2',[],n_pos)',2));
    shrink = max(1 - mu./max(rows_norm,mu),0);
    active_set = (shrink > 0);    
    shrink = shrink(active_set);
    if n_orient>1
        active_set = repmat(active_set,1,n_orient);
        active_set = reshape(active_set',n_orient,[]);
        active_set = active_set(:)';
    end
    temp = reshape(repmat(shrink,1,n_orient),length(shrink),n_orient)';
    temp = temp(:);
    Y = Y(find(active_set),:).*repmat(temp,1,size(Y,2));
end


function [Y active_set] = prox_l1(Y,lambda,n_orient)
    n_pos = size(Y,1)/n_orient;
    norms = sqrt(sum(reshape((abs(Y).^2),n_orient,[]),1));
    shrink = max(1-lambda./max(norms,lambda),0);
    shrink = reshape(shrink',n_pos,[]);
    active_set = logical(sum(shrink,2));    
    shrink = shrink(find(active_set),:);
    
    if n_orient>1
        active_set = repmat(active_set,1,n_orient);
        active_set = reshape(active_set',n_orient,[]);
        active_set = active_set(:)';
    end
    
    Y = Y(active_set,:);
    if length(Y) > 0
        for i = 1:n_orient
            Y(i:n_orient:size(Y,1),:) = Y(i:n_orient:size(Y,1),:).*shrink;
        end
    end
end

% function a = norm_l21(Z)
%     
%     if isempty
% end
% function a = norm_l1(Z)
% 
% end

function k = lipschitz_contant(y, L, tol, a, M)
    Nt = size(y,2);
    Nd = size(L,2);
    iv = ones(Nd,Nt);
    v = dgtreal(iv', 'gauss', a,M);
    
    T = size(v,2);
    K = size(v,1);
    
    l = 1e10;
    l_old = 0;
    fprintf('Lipschitz constant estimation: \n')
    rev_line = '';
    for i = 1 : 100               
       msg = sprintf('Iteration = %d, Diff: %d, Lipschitz Constant: %d ',i,abs(l-l_old)/l_old,l);
       fprintf([rev_line, msg]);
       rev_line = repmat(sprintf('\b'),1,length(msg));
       l_old = l;
       aux = idgtreal(v,'gauss',a,M)';
       iv = real(aux);
       Lv = L*iv;
       LtLv = L'*Lv;
       w = dgtreal(LtLv', 'gauss', a,M);
       l = max(max(max(abs(w))));
       v = w/l;
       if abs(l-l_old)/l_old < tol
           break
       end
    end
    fprintf('\n');
    k = l;

end