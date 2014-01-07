function y = nip_addnoise(L,J,t,snr)
% y = nip_addnoise(L,J,t,snr)
%   Input:
%       L -> NcxNd. Leadfield matrix
%       J -> NdxNt. Clean brain activity
%       snr -> scalar. Resulting SNR
%   Output:
%       y -> NcxNt. Noisy EEG signal
% 
%
% Juan S. Castano C.
% 2 Jan 2013
Nd = size(L,2)/3;
Nt = length(t);
Nc = size(L,1);

% Central frequency for "noise wavelets"
interval = [1 30];
freq = interval(1) + (interval(2)-interval(1)).*rand(round(Nd*0.15),1);

% Time shift for "noise wavelets"
interval = [t(1) t(end)];
p_shift = interval(1) + (interval(2)-interval(1)).*rand(round(Nd*0.15),1);

% Active sources
Nnoise = length(freq);
idx = randsample(Nd, Nnoise);
dir = randn([Nd, 3]);
dir = dir.*repmat(sqrt(sum(dir.^2,2)),1,3);

Jnoise = zeros(3*Nd,Nt);

n = 1;
for i = idx'
    idxr = ((i-1)*3+1):((i-1)*3+3);
    series = aux_wavelet(freq(n),t,p_shift(n));
    Jnoise(idxr,:) = dir(n,:)'*series;
    n = n+1;
end

brainnoise = L*Jnoise;
brainnoise = brainnoise/norm(brainnoise,'fro');

sensornoise = randn(Nc,Nt);
sensornoise = sensornoise/norm(sensornoise,'fro');

beta = 0.6;
noise = beta*brainnoise + (1-beta)*sensornoise;


signal = L*J;
signal = signal/norm(signal,'fro');

noise = noise/norm(noise,'fro');
 
snr = 10^(snr/10);
y = signal + (1/snr)*noise;
y = y*norm(L*(J),'fro');

end

function series = aux_wavelet(f0,t,phase_shift)



% Normalization terms for the wavelet
sigma_f = f0/7;
sigma_t = 1/(2*pi*sigma_f);

% "source" contains the time courses of active sources
series =  real(exp(2*1i*pi*f0*t).*...
    exp((-(t-phase_shift).^2)/(2*sigma_t^2)));
end