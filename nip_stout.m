function [J_rec extras] = nip_stout(y,L,B,varargin)
% function [J_rec extras] = nip_stout(y,L,B,options)
% Input:
%         y -> NcxNt. Matrix containing the data,
%         L -> Ncx3Nd. Lead Field matrix
%         B -> NdxNk. Spatial basis functions
%         options -> struct.
%                 options.spatial_reg -> scalar. Percentage of spatial
%                       regularization (between 0 and 1). 
%                 options.temp_reg -> scalar. Percentage of temporal
%                       regularization (between 0 and 1).
%                 options.a -> scalar. Time shift for the time frequency
%                       transform.
%                 options.m -> scalar. Frequency bins for the time frequency
%                       transform.
%                 options.tol -> Scalar. Default 1e-3
%   Output:
%         J_rec -> 3NdxNt. Reconstructed activity (solution)
%         extras. -> Currently empty
%
% Juan S. Castano C.
% jscastanoc@gmail.com
% 19 Aug 2013

Ltemp = nip_translf(L);
for i = 1:3
    Lnew(:,:,i) = Ltemp(:,:,i)*B;
end
Lnew = nip_translf(Lnew);
[J_r,~] = nip_tfmxne_port(y,Lnew,varargin{:});
J_rect = nip_translf(J_r');
J_rect = permute(J_rect,[2 1 3]);
for i = 1:3
    J_rec(:,:,i) = B*J_rect(:,:,i);
end
J_rec = permute(J_rec,[2 1 3]);
J_rec = nip_translf(J_rec)';
extras =[];
end