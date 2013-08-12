function [J_est, extras] = nip_sflex(y, L, basis, reg_par)
%  [J_est, extras] = nip_sflex(y, L, basis)
% Implements "Large-scale EEG/MEG source localization with spatial
% flexibility." by Haufe et al 2011
%
% Input:
%       y -> NcxNt. Matrix containing the data,
%       L -> NcxNd. Lead Field matrix
%       basis-> NdxNs. Matrix containing the spatial basis functions.
%       reg_par-> Scalar. Regularization parameter (1e-6 by default)
% Output:
%       J_rec -> NdxNt. Reconstructed activity (solution)
%       extras.regpar -> Scalar. Optimum regularization parameter
%
% Additional comments: Uses the DAL optimization toolbox.
% 
% Juan S. Castano C.
% 13 June 2013


[Nc Nt] = size(y);
Nd = size(L,2);

if nargin <=3
    [~,extras]=nip_loreta(y,L,speye(Nd));
    reg_par = extras.regpar*100;
end
% reg_par = 1e-2;

% --- Not sure if I should do this here (you'll need a lot of ram if data
% is big enough
A = sparse(kron(speye(Nt), L*basis)); 

nbasis = size(basis,2);


[xx,status]=dalsqgl(zeros(nbasis/3,3*Nt), A, y(:), reg_par);
if sum(xx) == 0
    [~,extras]=nip_loreta(y,L,speye(Nd));
    reg_par = extras.regpar*30;
    [xx,status]=dalsqgl(zeros(nbasis/3,3*Nt), A, y(:), reg_par);
end
J_est = basis*reshape(xx,nbasis,Nt);
extras =[];

end