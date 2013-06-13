function [y_proj, y_rec, Ur, Er] = nip_tempcomp(y, t, pf, ev)
%  [y_proj, y_rec, Ur, Er] = nip_tempcomp(y, t, pf, ev)
% This function projects (and also reconstructs) data y (EEG) using the
% main components of its covariance matrix of the DCT of the data (yy'). It can
% also be used to filter the data using the pf parameters which are the
% cut-off frequencies of a band-pass filter.
%
% Input:
%       y -> NcxNt. Matrix containing the data,
%       t -> 1xNt. Time vector (seconds)
%       pf -> 1x2. Vector with low cut-off frequency and high cut-off frequency (default 0 60).
%       ev -> scalar. If < 1, returns uses the necessary number of
%       eigenvectors to achieve ev*100 % of explained variance. If > 1
%       indicates the number of eigenvectors to use.
%
% Output:
%       y_proj -> NcxNr. Project data (Nr is the number of components used)
%       y_rec  -> NcxNt. Reconstructed data
%       Ur     -> NtxNr. Projection matrix (Eigenvectors in the time domain).
%       Er     -> 1xNr.  Eigenvalue associated with each eigenvector in Ur
%
% Additional comments:
%       Parts of this function are based on a script written by:
%               Jose David Lopez - ralph82co@gmail.com
%				Gareth Barnes - g.barnes@fil.ion.ucl.ac.uk
%				Vladimir Litvak - litvak.vladimir@gmail.com
%
%
% Juan S. Castano C.
% 7 Jun 2013




[Nc, Nt] = size(y);
T = dctmtx(Nt); % Transformation matrix for the DCT
fs = 1/(t(2)-t(1)); % Sample frequency

dct = t*fs+1;
dct = dct/2/t(end);

% Filter in the frequency domain
if nargin >=3 && ~isempty(pf)
    lpf = pf(1);
    hpf = pf(2);    
    j = find( (dct >= lpf) & (dct <= hpf) );
else
    j = 1:Nt;
end

T = T(:,j);
dct = dct(j);

yty = T'*y'*y*T;

[U E] = svd(yty);

if nargin <4
    ev = 10;
end

if ev > 1
    Nr = ev;
else
    cum_var = cumsum(diag(E))/sum(diag(E));
    idx = find(cum_var >= ev);
    Nr =min(idx); % Number of temporal modes
end

Ur = U(:,1:Nr);
Er = diag(E(1:Nr,1:Nr));

pc_time = T*Ur; %Proyection matrix


y_proj = y*pc_time;

y_rec = y_proj*pc_time';

Ur = pc_time;
end