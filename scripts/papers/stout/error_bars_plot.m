% Error bars

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 5 sep  2013           %
%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization and creation of the structures used in some functions
close all; clc; clear

Ntrials = [5 20 250];
act_sources = [5];
% Ntrials = [210];
snr_bio = 5;
% snr_bio = [5]
snr = [];
% t = model.t;
% L = model.L;
% fs = model.fs
% cortex = model.cortex;
% clear model;

n_exp = 50;


jobs_c = 1;
% methods = {'TF-MxNE','S+T'};
methods = {'TF-MxNE','S+T','LOR','S-FLEX'};
% methods = {'S-FLEX'};
errors = [1,2,3,4,5];
dummy_counter = 0;
% total_iter = length(errors)*numel(methods)*(n_exp)*length(act_sources)*length(Ntrials);
lspec = {'b','k','g','r'};
% ff = figure;

for i = 1:5
    ff(i) = figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);;
end
for c_meth = 1:numel(methods)
    nn = 1;
    for i = Ntrials
        for j = 1:n_exp
            cur_jobs = [];
            copy_res = {};
            for l = 1:length(act_sources)
                dir = strcat('/mnt/data/error/montreal_sampleall_false/',num2str(act_sources(l)));
                file_name = strcat(dir,'/',methods{c_meth},'Exp',num2str(j),'Ntrials',...
                    num2str(i),'BioNoise',num2str(snr_bio),'.mat');
                load(file_name,'er');
                for ier = 1:5
                    error_mat{ier}(j,nn) = er(ier);
                end
            end
        end
        nn = nn+1;
    end
    
    for ier = 1:5
        inter_vec = (Ntrials(1):5:Ntrials(end));
        figure(ff(ier))
        aux = sort(error_mat{ier},1,'descend');
        %         plot(inter_vec,interp1(Ntrials,mean(aux(10:end,:)),inter_vec,'cubic'),lspec{c_meth})
        semilogx(inter_vec,interp1(Ntrials,mean(aux),inter_vec,'cubic'),lspec{c_meth})
        xlim([inter_vec(1) inter_vec(end)])
        hold on
    end
end
err = {'emd','corr','geo','wgeo','ev'};
titles = {'EMD', 'Unnormalized correlation a.u.', 'Geodesic Distance' , 'Weighted Geodesic distance' , 'Explained Variance (%)'};
for ier = 1:5
    figure(ff(ier))
    if ier == 1
        legend(methods,'Location','East')
    else
    legend(methods,'Location','SouthEast')
    end
    xlabel('Ntrials')
    ylabel(titles{ier})
    errorname = strcat('ActS',num2str(act_sources),'_',err{ier});
    fig_name = strcat('/mnt/data/ResultadosBerlin/FiguresOthers/',errorname);   
    savefig(fig_name, ff(ier), 'eps');
end