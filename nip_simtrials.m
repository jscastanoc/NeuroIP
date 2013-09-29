function [y, Jclean, actidx] = nip_simtrials(L, dip_pos, CleanBrain, t, Nsp, Ntrials,snr_meas, snr_bio, options)
% [y, Jclean, actidx] = nip_simtrials(L, dip_pos, act, t, Nsp, Ntrials,snr_meas, snr_bio)
%  Input:
%         L -> Ncx3Nd. Lead Field matrix
%         dip_pos -> Ndx3. Coordinates of each dipole
%         CleanBrain -> NdxNt. Brain activity without noise
%         t -> 1xNt. Time vector (in secs).
%         Nsp -> scalar. Number of Spurious/noise sources.
%         Ntrials -> scalar. Number of trials to simulate.
%         snr_meas -> scalar. SNR at the sensor level in each trial.
%         snr_bio -> scalar. SNR at the source level in each trial.
%         options -> same as nip_simulate_activity
%   Output:
%         y -> NcxNt. Measurements averaged across trials.
%         Jclean -> 3NdxNt. Ground true. Brain activity without noise.
%         actix. -> Nactx1. Indices of the active dipoles (as indexed in
%               dip_pos.
%
% Juan S. Castano C.
% jscastanoc@gmail.com
% 19 Aug 2013


rng('default')
rng('shuffle')
[Nc Nd] = size(L);
Nt = length(t);


% [Jclean, actidx] = nip_simulate_activity(dip_pos, Nact, act, dir, t,options);
Jclean = CleanBrain.Jclean;
actidx = CleanBrain.actidx;
Nact = length(actidx);
mask = ones(size(L,2),Nt);
for i =1:Nact
   mask(actidx(i):actidx(i)+2,:) = zeros(3,Nt);
end

J = zeros(Nd,Nt);
y_avg = zeros(Nc,Nt);
fprintf('Simulating trials and averaging ...\n')
rev_line = '';

dir = randn(Nsp,3);
sp_act = mkpinknoise(Nt,Nsp)';
options.sample_all = true;
J = nip_simulate_activity(dip_pos,Nsp, sp_act, dir, t,options);
for i = 1:Nact
    J(actidx(i):actidx(i)+2,:) = zeros(3,Nt);
end

sp_scale = norm(full(Jclean))/(10^(snr_bio/20)*norm(J));
for n = 1:Ntrials
    if mod(n,10) == 0   
        msg = sprintf('Trial # %d',n);
        fprintf([rev_line, msg]);
        rev_line = repmat(sprintf('\b'),1,length(msg));
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation of the biological noise %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Number of "spurious" active sources
    sp_act = mkpinknoise(Nt,Nsp)';
    
    options.sample_all = true;
    dir = randn(Nsp,3);
    J = nip_simulate_activity(dip_pos,Nsp, sp_act, dir, t,options);
    for i = 1:Nact
        J(actidx(i):actidx(i)+2,:) = zeros(3,Nt);
    end
    J = sp_scale*J + Jclean;
    
    y = L*J;
    y = nip_addnoise(y,snr_meas);    
    y_avg = y_avg + y/Ntrials;
    
end
fprintf('\n')

y = y_avg;
Jclean = sparse(Jclean);
end
