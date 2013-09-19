function [J_est, extras] = nip_sflex(y, L, basis, reg_par)
%  [J_est, extras] = nip_sflex(y, L, basis, reg_par)
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

index = (1:3:Nd);
for i = 0:2
    L(:,index+i) = L(:,index+i)*basis; % J simulado FINAL
end
% --- Not sure if I should do this here (you'll need a lot of ram)
% L = nip_translf(L);
A = sparse(kron(speye(Nt), L)); 
nbasis = size(basis,2);

reg_par = 1e-4;
[xx,status]=dalsqgl(zeros(nbasis*3,Nt), A, y(:), reg_par);
for i = 0:2
    xx(index+i,:) = basis*xx(index+i,:); % J simulado FINAL
end
J_est = xx;

extras =[];

end