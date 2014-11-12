function h = nip_kalman_hyper(y,L, varargin)
% h = nip_kalman_hyper(y,L, varargin)
% Input:
%       y -> NcxNt. Matrix with Nc channels and Nt time points
%       L -> NdxNt. Original lead field matrix multiplied by the spatial
%       dictionaries (L = L*D)
%
%
% Output:
%       h -> NkxNt. Time series of the hyperparameters 
% TODO:
% Setup varargin
% Juan S. Castano
% jscastanoc@gmail.com
% 7 Oct 2013


% Measurement matrix for covariance matrix -> cov(y) = kron(L,L)Q(h)
[Nc,Nk] = size(L);
Nt = size(y,2);
A = zeros(Nc^2,Nk);
for i = 1:Nk
    auxL = kron(L(:,i),L(:,i));
    A(:,i) = auxL;
end
L = A;
clear A;
clear auxL;

% Compute measurement covariance matrices of data  with overlapped windows
window = -floor(Nc/24):1:floor(Nc/24);
% window = -floor(50):1:floor(50);
ycov = zeros(Nc^2,Nt);
for i = 1:Nt
    w_act = i + window;
    idx_del = find((w_act < 1) | (w_act >Nt));
    w_act(idx_del) = [];
    aux = y(:,w_act)*y(:,w_act)';
    ycov(:,i) = aux(:);
end

h = kalman_priors(ycov,L);
end

function h = kalman_priors(ycov,L)
[Nc,Nt] = size(ycov);
[Nk] = size(L,2);
h = zeros(Nk,Nt);
P = 1*speye(Nk);
A = 0.9*speye(Nk);
Q = 1*speye(Nk);
S = 1e4*speye(Nc);
Sinv = inv(S);
eyeNk = speye(Nk);

fprintf('Computing kalman... \n');
rev_line = '';
eta = 0;
etat = eta;
for i = 2:Nt
    tic
    msg = sprintf('Iteration # %u of %u\nElapsed time for this iteration: %f \nTotal time: %f \n',i,Nt,eta,etat);
    fprintf([rev_line, msg]);
    rev_line = repmat(sprintf('\b'),1,length(msg));
    hap = A*h(:,i-1);
    Pap = A*P*A'+Q;
    y_res = ycov(:,i) - L*hap;
    %    S_res = L*Pap*L' + S;
    S_resInv = Sinv - (Sinv*L/(inv(Pap+1e-5)+L'*Sinv*L))*L'*Sinv;
    K = Pap*L'*S_resInv;
    h(:,i) = hap + K*y_res;
    if sum(find(isnan(h)))
        error('Kalman Filter did not converge\nFailed at time instant %i of %i \n',i,Nt)
        return
    end
    eta = toc;
    etat = etat + eta;
    P = (eyeNk - K*L)*Pap;
end
end