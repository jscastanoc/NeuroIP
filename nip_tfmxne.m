function [J_rec extras] = nip_tfmxne(y,L,options)
% function [J_rec extras] = nip_tfmxne()
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


spatial_reg = options.spatial_reg;
temp_reg = options.temp_reg;

a = options.a;
M = options.m;  


c = dgtreal(y','gauss',a,M);
T = size(c,2);
K = size(c,1);
Z = sparse(Nd,K*T);
Y = Z;
J_rec = sparse(Nd,Nt);
J_recf = sparse(Nd,M);


tau = 1;
% lambda = 0.001; % Time regularization parameter
base_par= 0; % Spatial regularization parameter
tempGY = L'*y;
for i = 1:Nd
    basepar = max(norm(tempGY(i,:),2),base_par);
end
mu = basepar*spatial_reg; % Spatial regularization parameter
lambda = basepar*temp_reg; % Time regularization parameter
clear tempGY;

fprintf('Running TF-MxNE algorithm... \n');
for i = 1:5
    tic
   fprintf('Iteration # %d ',i)
   % line 5 of algorithm 1 (see Ref paper)
   Z_0 = Z;    
   
   % line 6 of algorithm 1 (see Ref paper)
   [act_dip,~] = ind2sub(size(Y),find(Y));
   if ~isempty(act_dip)
       for j = unique(act_dip)'
          aux = idgtreal(reshape(full(Y(j,:)'),K,T),'gauss',a,M);
          J_rec(j,:) = aux(1:Nt);
       end
   end
   error = y - L*J_rec;
   err_trans = dgtreal(error','gauss',a,M);
   err_trans = permute(err_trans,[3 1 2]);  
   err_trans = reshape(err_trans,Nc,[]);   
   arg_prox = (Y + mu*L'*err_trans);
   clear err_trans
   Z = prox(arg_prox,mu,lambda);
   clear arg_prox
   
   % line 7 of algorithm 1 (see Ref paper)
   tau_0 = tau;
   
   % line 8 of algorithm 1 (see Ref paper)
   tau = (1+sqrt(1+4*tau^2))/2;
   
   % line 9 of algorithm 1 (see Ref paper)
   Y = Z + ((tau_0-1)/tau)*(Z-Z_0);
   
   figure
   imagesc(full(abs(Y)))
   pause(0.01);
   fprintf('Elapsed time: %f seconds\n',toc)
end
[act_dip,~] = ind2sub(size(Z),find(Z));

fprintf('Transforming solution to the time domain: \n%d non-zero time series \n '...
    , length(unique(act_dip)))
for j = unique(act_dip)'
    J_recf(j,:) = idgtreal(reshape(full(Z(j,:)),K,T),'gauss',a,M);
end
J_rec = J_recf(:,1:Nt);
extras = [];
end

function Z = prox(Y,mu,lambda)

[P,K] = size(Y);

zpk = sparse(P,K);

% index = 1:3:P;
% n = 1;
% for p = index
%         normY(n,:) = sqrt(sum(Y(p:p+2,:).^2));
%         n = n+1;
% end
% normY = reshape(repmat(normY(:),3,1),P,K);
% A = Y./normY;

A = Y./abs(Y);
A(find(isnan(A))) = 0;
B = max(abs(Y)-lambda,zpk);
AB = A.*B;
clear A B;


% parfor p = 1:P
%     dip_num = ceil(p/3)
%     C(p) = sqrt(sum(max((sqrt(norm(Y(dip_num,:))^2 + ...
%         norm(Y(dip_num+1,:))^2 + norm(Y(dip_num+2,:))^2)-lambda).^2,0),2));
%     aux = mu./C(p);
%     aux(find(isnan(aux))) = 0;
%     CC(p) = max(1 - aux,0);
% end

parfor p = 1:P
    C(p) = sqrt(sum(max( (abs(Y(p,:))-lambda).^2,0),2));
    aux = mu./C(p);
    aux(find(isnan(aux))) = 0;
    CC(p) = max(1 - aux,0);
end

Z = AB.*repmat(CC',[1,K]);

end

