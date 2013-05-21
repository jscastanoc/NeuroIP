function tai = nip_error_tai(y,L,Jrec)
% tai = nip_error_tai(y,L,Jrec)
%
% Computes an temporal accuracy index to evaluate the temporal quality of an
% inverse solution.
%
% Input:
%       y -> NcxNt. Matrix containing the data,
%       L -> NcxNd. Lead Field matrix
%       J_rec -> NdxNt. Reconstructed brain activity.
% Output:
%       tai -> scalar. Temporal Accuracy Index
%
% Additional Comments: See Belardinelli et al 2012 for further information about
% this error measurement (Source Reconstruction Accuracy of MEG and EEG...)
%
% Juan S. Castano C.
% 21 May 2013.

SSR  = sum(var((y - L*Jrec),0,2)); % Residual Variance
SST  = sum(var(y,0,2)); % Total variance
tai  = (SST - SSR)/SST;  

end