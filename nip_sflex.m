function [J_est, extras] = nip_sflex(y, L, basis, reg_par)
%  [J_est, extras] = nip_sflex(y, L, basis)
%
% Juan S. Castano C.
% 13 June 2013


[Nc Nt] = size(y);
Nd = size(L,2);

if nargin <=3
    reg_par = 1e-6;
end


% --- Not sure if I should do this here (you'll need a lot of ram if data
% is big enough
A = sparse(kron(speye(Nt), L*basis)); 

nbasis = size(basis,2);

[xx,status]=dalsql1(zeros(nbasis*Nt,1), A, y(:), reg_par);

J_est = basis*reshape(xx,nbasis,Nt);
extras =[];
end