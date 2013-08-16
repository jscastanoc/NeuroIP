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
%   Output:
%         J_rec -> 3NdxNt. Reconstructed activity (solution)
%         extras. -> Currently empty
%
% Comments: 
% Almost a literal translation of the function tf_mixed_norm_solver found in:
% https://github.com/mne-tools/mne-python/blob/0ac3ac1a1634673da013109f38daf4f162cae117/mne/inverse_sparse/mxne_optim.py
% Juan S. Castano C.
% jscastanoc@gmail.com
% 14 Aug 2013
%TODO:
% Convergence criteria
% lipschitz constant
 
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

fprintf('Running TF-MxNE algorithm... \n');
active_set = logical(sparse(1,Nd));
Y_time_as = [];
Y_as = [];
lipschitz_k = 100; % FIX THIS!!!;
mu_lc = mu/lipschitz_k;
lambda_lc = lambda/lipschitz_k;
for i = 1:options.iter
    tic
   fprintf('Iteration # %d ',i)
   % line 5 of algorithm 1 (see Ref paper)
   Z_0 = Z;   active_set0 = active_set;
   
   if (sum(active_set) < size(R,1)) && ~isempty(Y_time_as)
       GTR = L'*R/lipschitz_k;
       A = GTR;
       A(Y_as,:) = A(Y_as,:) + Y_time_as; 
       [~, active_set_l21] = prox_l21(A,mu_lc,3);
       idx_actsetl21 = find(active_set_l21);
       
       aux = dgtreal(GTR(idx_actsetl21,:)','gauss',a,M);
       aux = permute(aux,[3 1 2]);  
       aux = reshape(err_trans,sum(active_set_l21),[]);   
       
       B = Y(idx_actsetl21,:) + aux
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
   [Z, active_set_l21] = prox_l21(Z,mu_lc,3); %??   
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

    %    figure
    %    imagesc(full(abs(Y)))
    %    pause(0.01);

       Y_as = active_set0|active_set;

       Y_time_as = zeros(Nd,M);

       for j = find(Y_as)
           Y_time_as(j,:) = idgtreal(reshape(full(Y(j,:)),K,T),'gauss',a,M)';
       end
    %    Y_time_as2 = idgtreal(reshape(fliplr(full(Y(find(Y_as),:))),K,T,[]),'gauss',a,M);
       R = y - L(:, find(Y_as))*Y_time_as(find(Y_as),1:Nt);
       figure(10)
       plot(R')
       pause(0.01)
       fprintf('Elapsed time: %f seconds\n',toc)   
   end
end

fprintf('Done!... \n Transforming solution to the time domain: \n%d non-zero time series \n '...
    , length(active_set))
i=1;
for j = find(active_set)
       J_recf(j,:) = idgtreal(reshape(full(Z(i,:)),K,T),'gauss',a,M)';
       i = i+1;
end

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
