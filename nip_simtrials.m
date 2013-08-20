function [y, Jclean, actidx] = nip_simtrials(L, dip_pos, act, t, Nsp, Ntrials,snr_meas, snr_bio)

[Nc Nd] = size(L);
Nt = length(t);

Nact = size(act,1);
dir = randn(Nact,3);
options.sample_all = true;
[Jclean, actidx] = nip_simulate_activity(dip_pos, Nact, act, dir, t, options);


mask = ones(size(L,2),Nt);
for i =1:Nact
   mask(actidx(i):actidx(i)+2,:) = zeros(3,Nt);
end

J = zeros(Nd,Nt);
y_avg = zeros(Nc,Nt);
fprintf('Simulating trials and averaging ... \n')
rev_line = '';

dir = randn(Nsp,3);
sp_act = mkpinknoise(Nt,Nsp)';
J = nip_simulate_activity(dip_pos,Nsp, sp_act, dir, t,options);
for i = 1:Nact
    J(actidx(i):actidx(i)+2,:) = zeros(3,Nt);
end

sp_scale = norm(Jclean)/(10^(snr_bio/20)*norm(J));
for n = 1:Ntrials
    if mod(n,10) == 0   
        msg = sprintf('\rTrial # %d',n);
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

end
