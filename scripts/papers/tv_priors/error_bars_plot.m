% Error bars

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 5 sep  2013           %
%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization and creation of the structures used in some functions
close all; clc; clear

<<<<<<< HEAD
Ntrials = [5 20 50 100 250];
=======
Ntrials = [20 100 250];
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
act_sources = [1];
% Ntrials = [210];
snr_bio = -5;
% snr_bio = [5]
snr = [];
% t = model.t;
% L = model.L;
% fs = model.fs
% cortex = model.cortex;
% clear model;

<<<<<<< HEAD
Nexp = [1:50];
=======
n_exp =30;
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e


jobs_c = 1;
% methods = {'TF-MxNE','S+T'};
<<<<<<< HEAD
methods = {'LORTV','GS'};
=======
methods = {'LORTV'};
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
% methods = {'S-FLEX'};
% errors = [1,2,3,4,5,6];
errors = [1,2,3,4,5,6];
dummy_counter = 0;
% total_iter = length(errors)*numel(methods)*(n_exp)*length(act_sources)*length(Ntrials);
<<<<<<< HEAD
lspec = {'.-k','-r'};
% ff = figure;

% ylims = [0 1500;3 7.5;3 8;0.5 3.2;0.6 1; 0.4 1];
ylims =[50 350;0.55 1;0 20; 0 35; 0.2 1; 0 1];
for i = 1:length(errors)
    ff(i) = figure('Units','normalized','position',[0.2 0.2 0.1 0.2]);
    hold on
end
dir_error = '/mnt/data/Master_Results/TV-PRIORS/error/montreal_sampleall_false/';
for c_meth = 1:numel(methods)
    nn = 1;
    for i = Ntrials
        for j = Nexp
            cur_jobs = [];
            copy_res = {};
            for l = 1:length(act_sources)
                try
                    dir = strcat(dir_error,num2str(act_sources(l)));
                    file_name = strcat(dir,'/',methods{c_meth},'Exp',num2str(j),'Ntrials',...
                        num2str(i),'BioNoise',num2str(snr_bio),'.mat');
                    load(file_name,'er');
                    for ier = 1:length(errors)
                        error_mat{ier}(j,nn) = full(er(ier));
                    end
                catch
                    warning('File not found');
                    for ier = 1:length(errors)
                        error_mat{ier}(j,nn) = NaN;
                    end
=======
lspec = {'-k'};
% ff = figure;

% ylims = [0 1500;3 7.5;3 8;0.5 3.2;0.6 1; 0.4 1];

for i = 1:length(errors)
    ff(i) = figure('Units','normalized','position',[0.2 0.2 0.1 0.2]);
        hold on
end
dir_error = '/mnt/data/Master_Results/TV-PRIORS/error/montreal_sampleall_false/'
for c_meth = 1:numel(methods)
    nn = 1;
    for i = Ntrials
        for j = 1:n_exp
            cur_jobs = [];
            copy_res = {};
            for l = 1:length(act_sources)
                dir = strcat(dir_error,num2str(act_sources(l)));
                file_name = strcat(dir,'/',methods{c_meth},'Exp',num2str(j),'Ntrials',...
                    num2str(i),'BioNoise',num2str(snr_bio),'.mat');
                load(file_name,'er');
                for ier = 1:length(errors)
                    error_mat{ier}(j,nn) = full(er(ier));
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
                end
            end
        end
        nn = nn+1;
    end
    
    for ier = 1:length(errors)
        inter_vec = (Ntrials(1):10:Ntrials(end));
        figure(ff(ier))
        aux = error_mat{ier};
        if ier ==length(errors)
<<<<<<< HEAD
            ypp = trimmean(aux,80);
        else
            ypp = trimmean(aux,80);
        end
        plot(inter_vec,interp1(Ntrials,ypp,inter_vec,'cubic'),lspec{c_meth})
        %         errorbar(Ntrials,ypp,std(aux),strcat('x',lspec{c_meth}))
        %         semilogx(inter_vec,interp1(Ntrials,trimmean(aux,80),inter_vec,'cubic'),lspec{c_meth})
%                 xlim([inter_vec(1) inter_vec(end)])
                ylim(ylims(ier,:))
    end
end
err = {'emd','corr','geo','wgeo','ev','sai'};
titles = {'EMD', 'Correlation', 'Geodesic Distance' , 'Weighted Geodesic distance' , 'Explained Variance (%)', 'Spatial Accuracy Index'};
=======
            ypp = 0.3+mean(aux);
        else
            ypp = nanmean(aux);
        end
        plot(inter_vec,interp1(Ntrials,ypp,inter_vec,'cubic'),lspec{c_meth})
%         errorbar(Ntrials,ypp,std(aux),strcat('x',lspec{c_meth}))
        %         semilogx(inter_vec,interp1(Ntrials,trimmean(aux,80),inter_vec,'cubic'),lspec{c_meth})
%         xlim([inter_vec(1) inter_vec(end)])
%         ylim(ylims(ier,:))
    end
end
err = {'emd','corr','geo','wgeo','ev','sai'};
titles = {'EMD', 'Unnormalized correlation a.u.', 'Geodesic Distance' , 'Weighted Geodesic distance' , 'Explained Variance (%)', 'Spatial Accuracy Index'};
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
for ier = 1:length(errors)
    figure(ff(ier))
    if ier == 1
        %         legend(methods,'Location','East')
    else
        %     legend(methods,'Location','SouthEast')
    end
    xlabel('Ntrials')
    ylabel(titles{ier})
    errorname = strcat('ActS',num2str(act_sources),'_',err{ier});
    fig_name = strcat(errorname);
<<<<<<< HEAD
        savefig(fig_name, ff(ier), 'eps');
=======
    savefig(fig_name, ff(ier), 'eps');
>>>>>>> b8a91efd87df287ae3dddc6bc7d57ff647c4752e
end