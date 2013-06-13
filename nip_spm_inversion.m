function [J_est, extras] = nip_spm_inversion(y,L,Q,Qe)
% function [J_est, extras] = nip_spm_inversion(y,L,Q,Qe)
% This is a wrapper for the Inverse solution found in SPM
% 
%
% Additional Comments: This function is based on a script written by
%           	Jose David Lopez - ralph82co@gmail.com
%				Gareth Barnes - g.barnes@fil.ion.ucl.ac.uk
%				Vladimir Litvak - litvak.vladimir@gmail.com
%
% Juan S. Castano
% 9 Jun 2013
YY = y*y';
Qml = {Qe, L*Q*L'};
V = exp(-2)*trace(YY)*eye(size(y,1))/size(y,1);

[Cy,h,Ph,F] = spm_reml_sc(YY,[],Qml,1,-4,16,V);

Cp  =  h(2)*diag(Q);
LCp =  h(2)*L*Q;

M     = LCp'/Cy;

% Q = h(2)*Q;
% aux = Q*L';
% M = aux/(L*aux+h(1)*Qe);

J_est = M*y;

Qpost = J_est*J_est';
Qml = {Qe, L*Qpost*L'};
[Cy,h,Ph,F] = spm_reml_sc(YY,[],Qml,1,-4,16,V);

extras.par = h;
extras.FE = F;