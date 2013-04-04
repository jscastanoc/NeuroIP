function y = nip_addnoise(x,snr)
% x = nip_addnoise(x,snr)
%   Input:
%       x -> NcxNt. Channel signals
%       snr -> scalar. Resulting SNR
%   Output:
%       y -> NcxNt. Signals with added white noise
% 
% Additional Comments:
% Based on a script written by Prof. Gareth Barnes
%
% Juan S. Castano C.
% 1 Feb 2013

allchanstd=(std(x'));
meanrmssignal=mean(allchanstd);
for i = 1:size(x,1)
    ch_noise = meanrmssignal.*randn(size(x(i,:)))/(10^(snr/20));
    allchannoise(i,:)=ch_noise;
    y(i,:) = x(i,:) + ch_noise;
end