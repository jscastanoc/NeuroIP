close
cfg.fs = 100; % Sample frequency (Hz)
t = 0:1/cfg.fs:1.5; % Time vector (seconds)

j = randi([1,50],1,1);
k = randsample([5 20 50 100 250],1);
i = randsample([-15 -5],1);

Nf = 3,
dir = strcat('D:/Datasets/sim_trials_final/',num2str(Nf),'/');
file_name = strcat(dir,'Exp',num2str(j),'Ntrials',...
    num2str(k),'BioNoise',num2str(i),'.mat');
load(file_name);
y = mgjob.results{1};
Jclean = mgjob.results{2};
actidx = mgjob.results{3};


% idx = [];
act_fig = figure('Units','normalized','position',[0.1 0.1 0.2 0.3]);
lspec = {'-k','*-k','x-k','d-k','o-k'};
% for i = 1:length(actidx)    
    plot(t,Jclean(actidx*3,:))
    hold on
%     idx(:,i) = ((actidx(i)-1)*3+1:(actidx(i)-1)*3+3);
% end
ylabel('Amplitude (a.u.)')
xlabel('Time (sec)')
ylim([-1 1])
savefig(strcat('ActN',num2str(Nf)),act_fig,'eps')
% idx = idx(:);
% for i 