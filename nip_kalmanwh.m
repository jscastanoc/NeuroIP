function [x, extras] = nip_kalmanwh(y, L, Laplacian,par)
% [x, extras] = nip_kalmanwh(y, L, Laplacian)
% Computes the kalman filter estimation using a whitening to eliminate
% spatial correlations (Presented in Galka et al. 2004). See Documentation
% for further details.
%
% Input: 
%       y -> NcxNt. Matrix containing the data,
%       L -> NcxNd. Lead Field matrix
%		Laplacian -> NdxNd. Spatial Laplacian matrix (information about dipole Neighbors.
%       par -> 5x1 (optional). Vector containing the parameters a1 b1 a2
%       (parameters for the dynamical model), and sigma and epsi (variances)
%
% Output:
% 		x -> NdxNt. Estimated neural activity.
%		extra -> Struct. Contains extra output parameters:
%               Currently empty
%
% Additional Comments:
%
% Juan S. Castano C.
% 17 Feb 2013

if nargin <=3
    aic_fun = @(par) aic(y,L,Laplacian,par);
    % options = optimset('Display','iter','tolX',1e-6);
    options = optimset('Display','iter','tolX',1e-4,'TolFun',1e-2);
    par = fminsearch(aic_fun,[0.9 0.001 0.001 1e-3 1e-4],options);
elseif nargin == 4
    par;
end
a = par(1:3);
sigma = par(4);
epsi = par(5);

Nt = size(y,2);
Nc = size(L,1);
Nd = size(L,2);

eye_2 = eye(2);
eye_Nc = speye(Nc);

K = L*inv(Laplacian);

Akf = [a(1) a(3); 1 0];
Bkf = [a(2) 0   ; 0 0];

Qv = zeros(Nc,2);
covSigmaMat = [sigma^2 0; 0 0];

[Jwth, ~] = nip_loreta(y(:,1:2),L, inv(Laplacian'*Laplacian));
Jwth = Laplacian*Jwth;
P = 2*ones(Nd,1);
Pwth_ap = zeros([2,2,Nd]);
for k = 3:Nt
    R = zeros([Nc,Nc]);    
    LJwth = Laplacian*Jwth(:,k-1);
    for i = 1:Nd        
        Qv(:,1) = K(:,i);
        Jwth_ap(:,i) = Akf*[Jwth(i,k-1); Jwth(i,k-2)] + ...
            Bkf*[LJwth(i);0];
        Pwth_ap(:,:,i) = Akf*P(i)*Akf'+ covSigmaMat;          
        R = R + Qv*Pwth_ap(:,:,i)*Qv';
    end    
    R = R + epsi^2*eye_Nc;
    Y_pre = K*Jwth_ap(1,:)';
    dY = y(:,k)-Y_pre;
    Rinv = inv(R);
    for i = 1:Nd        
        Qv(:,1) = K(:,i);
        G = Pwth_ap(:,:,i)*Qv'*Rinv;
        aux = Jwth_ap(:,i) + G*dY;
        Jwth(i,k) = aux(1);
        aux = (eye_2 - G*Qv)*Pwth_ap(:,:,i);
        P(i) = aux(1,1);
    end
end
x = inv(Laplacian)*Jwth;
% extras.P = Not implemented;
extras.par = par;

function val = aic(y, L, Laplacian,par)


a = par(1:3);
sigma = par(4);
epsi = par(5);

Nt = size(y,2);
Nc = size(L,1);
Nd = size(L,2);

eye_2 = eye(2);
eye_Nc = speye(Nc);

K = L*inv(Laplacian);

Akf = [a(1) a(3); 1 0];
Bkf = [a(2) 0   ; 0 0];

Qv = zeros(Nc,2);
covSigmaMat = [sigma^2 0; 0 0];
covObs = epsi^2*eye_Nc;

[Jwth, ~] = nip_loreta(y(:,1:2),L, inv(Laplacian'*Laplacian));
Jwth = Laplacian*Jwth;
P = 2*ones(Nd,1);
Pwth_ap = zeros([2,2,Nd]);
logL = 0;
for k = 3:Nt
    R = zeros([Nc,Nc]);    
    LJwth = Laplacian*Jwth(:,k-1);
    for i = 1:Nd        
        Qv(:,1) = K(:,i);
        Jwth_ap(:,i) = Akf*[Jwth(i,k-1); Jwth(i,k-2)] + ...
            Bkf*[LJwth(i);0];
        Pwth_ap(:,:,i) = Akf*P(i)*Akf'+ covSigmaMat;          
        R = R + Qv*Pwth_ap(:,:,i)*Qv';
    end    
    R = R + covObs;
    Y_pre = K*Jwth_ap(1,:)';
    dY = y(:,k)-Y_pre;
    Rinv = inv(R);
    for i = 1:Nd        
        Qv(:,1) = K(:,i);
        G = Pwth_ap(:,:,i)*Qv'*Rinv;
        aux = Jwth_ap(:,i) + G*dY;
        Jwth(i,k) = aux(1);
        aux = (eye_2 - G*Qv)*Pwth_ap(:,:,i);
        P(i) = aux(1,1);
    end    
    if k>Nt/2
        logL = logL + log(abs(det(R)))+dY'*Rinv*dY+Nc*log(2*pi);
    end
end
val = logL + 10;