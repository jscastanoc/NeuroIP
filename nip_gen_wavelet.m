function series = nip_gen_wavelet(f0,t,phase_shift)
% series = nip_gen_wavelet(f0,t,phase_shift)
% Generate real morley wavelet along vector time t with central frequency f0 and
% phase shift phase_shit
% Normalization terms for the wavelet
sigma_f = f0/7;
sigma_t = 1/(2*pi*sigma_f);

% "source" contains the time courses of active sources
series =  real(exp(2*1i*pi*f0*t).*...
    exp((-(t-phase_shift).^2)/(2*sigma_t^2)));
end