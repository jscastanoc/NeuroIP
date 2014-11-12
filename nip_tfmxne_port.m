function [J_rec extras] = nip_tfmxne_port(y,L,varargin)
% function [J_rec extras] = nip_tfmxne_port(y,L,varargin)
% Solve the inverse problem using the TF-MxNE approach as introduced by
% Gramfort et al 2012.
%  Input:
%         y -> NcxNt. Matrix containing the data,
%         L -> Ncx3Nd. Lead Field matrix
%         Additional options -> Key - Value pair:
%                 'sreg' -> scalar. Percentage of spatial
%                       regularization (between 0 and 1).
%                 'treg' -> scalar. Percentage of temporal
%                       regularization (between 0 and 1).
%                 a -> scalar. Time shift for the time frequency
%                       transform.
%                 m -> scalar. Frequency bins for the time frequency
%                       transform.
%                 tol -> Scalar. Default 5e-2
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

global log_file

fname = '/home/jscastanoc/Desktop/dev_logger_matlab.log';
log_file = fopen(fname,'w');
fprintf(log_file,'Starting fista \n')



Ndor = size(L,2);
idx = 1:1:Ndor;
L_or = L;

L = L(:,idx);
[Nc Nd] = size(L);
Nt = size(y,2);

p = inputParser;

def_a = 8;
def_m = 112;
def_sreg = 80;
def_treg= 1;
def_maxiter = 500;
def_tol = 1e-8;
def_resnorm = 0.3;
def_lipschitz = [];
def_optimres = false;
def_Winv = [];

addParamValue(p,'tstep',def_a);
addParamValue(p,'wsize',def_m);
addParamValue(p,'sreg',def_sreg);
addParamValue(p,'treg',def_treg);
addParamValue(p,'maxiter',def_maxiter);
addParamValue(p,'tol',def_tol);
addParamValue(p,'resnorm',def_resnorm);
addParamValue(p,'lipschitz',def_lipschitz);
addParamValue(p,'optimres',def_optimres);
addParamValue(p,'Winv',def_Winv)

parse(p,varargin{:})
options = p.Results;


tol = options.tol;
sreg = options.sreg;
treg = options.treg;
a = options.tstep;
M = options.wsize;
lipschitz_k = options.lipschitz;
maxiter = options.maxiter;




% Initialization of the TF-MxNE algorithm
c = dgtreal(y','sine',a,M);
T = size(c,2);
K = size(c,1);
Z = sparse(0,K*T);
Y = sparse(Nd,K*T);
J_rec = sparse(Nd,Nt);


% Calculate scale leadfield matrix to use normalized reg parameters
tau = 1;
tempGY = L'*y;
aux = sum(reshape(tempGY.^2',[],Nd/3)',2);
basepar = 0.01*sqrt(max(aux(:)))

L = L/basepar;
L_or = L_or;
clear tempGY;

R = y;

active_set = logical(sparse(1,Nd));
Y_time_as = [];
Y_as = [];

if isempty(options.lipschitz)
    lipschitz_k = 1.1*lipschitz_contant(y, L, 1e-3, a, M);
end
fprintf(log_file,'Lipschitz constant : %s \n', lipschitz_k)

mu_lc = sreg/lipschitz_k;
lambda_lc = treg/lipschitz_k;
stop =false;
fprintf('Running TF-MxNE algorithm... \n');

eta = 0;

rescum = [];

temp =  reshape(full(Z)',K,T,[]);
error = inf;
res_0 = inf;
nn = 1;
while true
    rev_line = '';
    Z = sparse(0,K*T);
    Y = sparse(Nd,K*T);
    J_rec = sparse(Nd,Nt);
    tau = 1;
    temp =  reshape(full(Z)',K,T,[]);
     active_set = logical(sparse(1,Nd));
    Y_time_as = [];
    Y_as = [];
    for i = 1:maxiter
        tic;        
        
        Z_0 = Z;   
        active_set0 = active_set;
        
        temp = dgtreal(R','gauss',a,M);
        temp = permute(temp,[3 1 2]);
        temp = reshape(temp,Nc,[]);
        
        Y = Y + L'*temp/lipschitz_k;
        [Z, active_set_l1] = prox_l1(Y,lambda_lc,3);
        [Z, active_set_l21] = prox_l21(Z,mu_lc,3);
        active_set = active_set_l1;
        active_set(find(active_set_l1)) = active_set_l21;
%         active_set = active_set_l1 & active_set_l21;
        
        error_0 = error;
        if norm(active_set - active_set0) == 0
            error = norm(abs(Z-Z_0))/norm(abs(Z_0));
        else
            error = inf;
        end

        
%         stop = ((sum(full(active_set)) > sum(full(active_set0))) && i > 1) ;
%         if stop
%             disp('Converged')
%             Z = Z_0;
%             active_set = active_set0;
%             break
%         end
        msg = sprintf('Iteration # %d, Stop: %d, Time per iter: %f \nDipoles~=0: %d \nRegPar S:%d T:%d \n'...
            ,i,error,eta,sum(full(active_set))/3, mu_lc*lipschitz_k, lambda_lc*lipschitz_k);
        fprintf([rev_line, msg]);
        rev_line = repmat(sprintf('\b'),1,length(msg));
        if i < options.maxiter 
            % line 7 of algorithm 1 (see Ref paper)
            tau_0 = tau;
            
            % line 8 of algorithm 1 (see Ref paper)
            tau = 0.5+sqrt(1+4*tau_0^2)/2;
            
            Y = sparse(Nd,K*T);
            dt = (tau_0-1)/tau;
            Y(find(active_set),:) = (1 + dt)*Z;
            Y(find(active_set0),:)= Y(find(active_set0),:) - dt*Z_0;
            
            
            Y_as = active_set0 | active_set;
            
%             temp =  reshape(full(Y)',K,T,[]);
            temp = full(Y(find(Y_as),:));
            temp =  reshape(temp',K,T,[]);
            temp = flipdim(temp,2);            
            temp = flipud(idgtreal(temp,'gauss',a,M))';
            Y_time_as = sparse(Nd,size(temp,2));
%             Y_old = Y_time_as;
            Y_time_as(find(Y_as),:) = temp;
            
            % Residual
            R = y - L(:, find(Y_as))*Y_time_as(find(Y_as),1:Nt);
        end
        
%         stop = error < tol || error_0 < error ||( sum(full(active_set))==0 && i > 1);
        stop = error < tol || ( sum(full(active_set))==0 && i > 1);
        if stop
            disp('Converged')
            break
        end
        eta = toc;
    end
    
    fprintf(' \nDone!... \nTransforming solution to the time domain: \n%d non-zero time series \n'...
        , sum(full(active_set)))
    
    temp =  reshape(full(Z)',K,T,[]);
    temp = flipdim(temp,2);
        
    temp = flipud(idgtreal(temp,'gauss',a,M))';
    J_recf = sparse(Nd,size(temp,2));
    J_recf(find(active_set),:) = temp;
    J_recf = J_recf(:,1:Nt);
    Jf = zeros(Ndor,Nt);
    Jf(idx,:) = J_recf;
    Jf = Jf*norm(y,2)/norm(L_or*Jf,2);
    
    
    resnorm = norm(y-L_or*Jf, 'fro')/norm(y, 'fro');
    fprintf('\nRESNORM = %8.5e\n',resnorm);
    
    if resnorm >=1
        J_rec = Jf_0;
    else
        Jf_0 = Jf;
        J_rec = Jf;
    end
%     if nn >= options.maxiter || resnorm < options.gof...
%             || ~options.optimgof || gof_0 < resnorm
    if nn >= options.maxiter || resnorm < options.resnorm...
            || ~options.optimres || resnorm >=1
        break;
    else
        mu_lc = 0.8*mu_lc;
        lambda_lc = 0.8*lambda_lc;
    end
    res_0 = resnorm;
    nn = nn+1;
end
if ~isempty(options.Winv)
    J_est = nip_translf(J_rec');
    clear J_rec
    for i = 1:3
        J_rec(:,:,i) = (Winv(:,:,i)*J_est(:,:,i)')';
    end
    J_rec= nip_translf(J_rec)';
end
extras.active_set = active_set;
end

% function sm =  safe_max_abs(A, ia)
% if isempty(ia)
%     sm = 0;
% else
%     sm = max(max(abs(A(ia))));
% end
% end
% 
% function sm = safe_max_abs_diff(A,ia,B,ib)
% if isempty(ia)
%     A = 0;
% else
%     A = A(ia);
% end
% if isempty(ia)
%     B = 0;
% else
%     B = B(ib);
% end
% if isempty(A)&&isempty(B)
%     sm = 0;
% else
%     sm = max(max(abs(A-B)));
% end
% end
% 

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
Y = Y(find(active_set),:);

if length(Y) > 0
    for i = 1:n_orient
        Y(i:n_orient:size(Y,1),:) = Y(i:n_orient:size(Y,1),:).*shrink;
    end
end
end

function k = lipschitz_contant(y, L, tol, a, M)

global log_file
% L = L(:,randsample(size(L,2),900));
Nt = size(y,2);
Nd = size(L,2);
iv = ones(Nd,Nt);
v = dgtreal(iv', 'sine', a,M);
T = size(v,2);
K = size(v,1);

l = 1e100;
l_old = 0;
fprintf('Lipschitz constant estimation: \n')
rev_line = '';
xx = tic;
for i = 1 : 100
    tic
    msg = sprintf('Iteration = %d, Diff: %d, Lipschitz Constant: %d\nTime per iteration %d ',i,abs(l-l_old)/l_old,l,toc);
    fprintf([rev_line, msg]);
    
    fprintf(log_file,'Lipschitz estimation: iteration = %d \n', i)
    rev_line = repmat(sprintf('\b'),1,length(msg));
    l_old = l;
    aux = idgtreal(v,'sine',a,M)';
    iv = real(aux);
    Lv = L*iv;
    LtLv = L'*Lv;
    w = dgtreal(LtLv', 'sine', a,M);
    l = max(max(max(abs(w))));
    v = w/l;
    if abs(l-l_old)/l_old < tol
        break
    end
end
fprintf('\n');
k = l;
toc(xx)
end