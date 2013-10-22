function h = nip_kalman_hyper_sparse(y,L, varargin)
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
window = -floor(Nc/2):1:floor(Nc/2);
% window = -floor(50):1:floor(50);
ycov = zeros(Nc^2,Nt);
for i = 1:Nt
    w_act = i + window;
    idx_del = find((w_act < 1) | (w_act >Nt));
    w_act(idx_del) = [];
    aux = y(:,w_act)*y(:,w_act)';
    ycov(:,i) = aux(:);
end

[h, lambda_est, time_tot] = RWL1_DF(ycov, matrix3d(L,Nt), matrix3d(eye(Nk),Nt));
end
