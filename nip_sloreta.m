function [J_est] = nip_sloreta(y,L)
% [J_est] = nip_sloreta(y,L)
% Calculate the inverse problem solution using the standarized LORETA
% solution
% Input:
%       y -> NcxNt. Matrix containing the data,
%       L -> NcxNd. Lead Field matrix
% Output:
%       J_rec -> NdxNt. Reconstructed activity (solution)
%       extras.regpar -> Scalar. Optimum regularization parameter
%
% Additional comments: See Grech et al. 2008 for further information.
%                   Note that with this method there isn't temporal info.
% Juan S. Castano C. 
% jscastanoc@gmail.com
% 15 May 2013


[Nc,Nd] = size(L);
Nt = size(y,2);
[Jlor, extras, T] = nip_loreta(y,L,eye(Nd));

S = T*L;
for i = 1:Nd
    J_est(i,:) = (1/S(i,i))*Jlor(i,:).^2;
end

end