function [x, extra] = nip_ira(y, L, Q, Laplacian, Np, covmat_obs,par,varargin)
% [x, extra] = nip_ira(y, L, Q, Laplacian, Np, covmat_obs,par)
% Computes the iterative regularization algorithm to estimate brain
% activity. See Documentation for further details.
% Input:
%       y -> NcxNt. Matrix containing the data,
%       L -> NcxNd. Lead Field matrix
%       Q -> NdxNd. Covariance matrix of the noise in the process equation (dipole level).
%		Laplacian -> NdxNd. Spatial Laplacian matrix (information about dipole Neighbors.
%		Np -> Scalar. Number of parameters for the dynamic model.
%			Options: 2 (linear first order model), 3 (linear second order model), 5 (nonlinear second order model).
%		covmat_obs -> NcxNc (optional). Covariance matrix of the noise in the observation equation (sensor level).
%       par      -> 2x1(Optional)
%
% Output:
% 		x -> NdxNt. Estimated neural activity.
%		extra -> Struct. Contains extra output parameters:
%               extra.mpar: Array containg the time series of the estimated
%               parameters.
% Additional Comments:
% 
% Juan S. Castano C.
% jscastanoc@gmail.com
% 27 Jan 2013

p = inputParser;
def_Winv = [];
addParamValue(p,'Winv',def_Winv);
parse(p,varargin{:})
options = p.Results;

if nargin < 6
   covmat_obs = [];
   par = [];
end

warning('This function is not fully tested, use it at your own risk')

Nd = size(L,2);
Nc = size(L,1);
Nt = size(y,2);

eye_Nd = speye(Nd);
eye_Nc = speye(Nc);

% WTW = inv(Q);
WTW = Q;

% Sensor noise covariance
if isempty(covmat_obs) 
	Sigma = eye_Nc;
else
	Sigma = covmat_obs;
end

% Parameter "noise" covariance
RTR = eye(Np);
RTR_inv = inv(RTR);


% Anonymous function for the process model and the derivative of such function with respect to parameters w
switch Np
    case 2
        f = @(xk_1, xk_2, w) (w(1)*eye_Nd + w(2)*Laplacian)*xk_1;
        G = @(xk_1, xk_2) [xk_1 Laplacian*xk_1];
    case 3
        f = @(xk_1, xk_2, w) (w(1)*eye_Nd + w(2)*Laplacian)*xk_1 +...
    (w(3)*xk_2);
        G = @(xk_1, xk_2) [xk_1 Laplacian*xk_1 xk_2];
    case 5
        f = @(xk_1, xk_2, w) (w(1)*xk_1  + w(2)*Laplacian*xk_1 +...
    w(3)*xk_1.^2 + w(4)*xk_1.^3 + w(5)*xk_2);
        G = @(xk_1, xk_2) [xk_1 Laplacian*xk_1 xk_1.^2 xk_1.^3 xk_2];
end

% Update model for the parameters w. The random term is to "estimulate" variation on the parameters.
% g = @(wk_1) (1-0.01*randn(1))*wk_1; 
% g = @(wk_1) 1*wk_1; 
g = @(wk_1) 1*wk_1; 


% Empirically set, this controls the "velocity" of change of the parameters along time
% gamma= 1;


% Optimize lambda

if ~isempty(par)
    lambda = par(1);
    gamma =par(2) ;
else
    gcv_fun = @(par) gcv(par,y,L,f,g,G,Q, RTR_inv, eye_Nc, eye_Nd, Sigma, Np);
    options = optimset('Display','iter','tolX',1e-5);
    par = fminsearch(gcv_fun, [0.001 0.01], options);
    lambda = par(1);
    gamma = par(2);
end

Lambda_k = (1/lambda)*Q;
Gamma_k = (1/gamma)*RTR_inv;

[x, ~] = nip_loreta(y(:,1:2),L, 'cov',Q,'Winv',options.Winv);
wk = 0*randn(Np,1);
% wk(1) = 0.9;
w_par=[];
fprintf('Computing IRA... \n');
rev_line = '';
eta = 0;
etat = eta;
for k = 3:Nt
    plot(w_par')
    pause(0.1)
    tic
    msg = sprintf('Iteration # %u of %u\nElapsed time for this iteration: %f \nTotal time: %f \n',k,Nt,eta,etat);
    fprintf([rev_line, msg]);
    rev_line = repmat(sprintf('\b'),1,length(msg));
    wk_1 = wk;
    %The parameters lambda and gamma could be updated online (in this
    %implementation we won't do that
    %     lambda(k) = lambda(k-1);
    %     gamma(k) = gamma(k-1);
    %     Lambda_k = lambda(k)*Q;
    %     Gamma_k = gamma(k)*RTR_inv;

    
    Gk = G(x(:,k-1), x(:,k-2));
    temp1 = Gk'*Lambda_k;
    wk = (temp1*Gk + Gamma_k)\(temp1*x(:,k-1) + Gamma_k*g(wk_1));
    
    temp1 = Lambda_k*L';
    x(:,k) = (eye_Nd - (temp1/(L*temp1 + Sigma))*L)*...
        ((temp1/Sigma)*y(:,k) + f(x(:,k-1),x(:,k-2), wk));
    
    w_par(:,k)=wk;
    
    eta = toc;
    etat = etat + eta;
end

if ~isempty(options.Winv)
    J_est = nip_translf(x');
    siJ = size(J_est);
    J_rec = permute(reshape(full(reshape(permute(J_est, [1 3 2]), siJ(1), [])*options.Winv), siJ(1), siJ(3), siJ(2)), [1 3 2]);
    J_rec = nip_translf(J_rec)';
end

x = J_rec*(norm(y,'fro')/norm(L*J_rec,'fro'));

extra.mpar = w_par;
extra.regpar = [lambda gamma];
end

% Generalized Cross Validation to obtain the regularization parameters
% (this function is minimized)
function gcv_val = gcv(par,y,L,f,g,G,Q, RTR_inv, eye_Nc, eye_Nd, Sigma, Np)

Nt = size(y,2);
Nd = size(L,2);
lambda = par(1);
gamma = par(2);

Gamma = (1/gamma)*RTR_inv;

Lambda = (1/lambda)*Q;
T = Lambda*L'/(L*Lambda*L' + Sigma);

den = trace(eye_Nc-L*T)^2;

x = zeros(Nd,Nt);
wk = zeros(Np,1);
for k = 3:Nt
    wk_1 = wk;
    
    Gk = G(x(:,k-1), x(:,k-2));
    temp1 = Gk'*Lambda;
    wk = (temp1*Gk + Gamma)\(temp1*x(:,k-1) + Gamma*g(wk_1));
    
    temp1 = Lambda*L';
    x(:,k) = (eye_Nd - (temp1/(L*temp1 + Sigma))*L)*...
        ((temp1/Sigma)*y(:,k) + f(x(:,k-1),x(:,k-2), wk));
end

% Take only the last half of the segment to allow transient to pass
num = norm(L*x(:,round(end/2):end)-y(:,round(end/2):end))^2;

gcv_val = num/den;

end
