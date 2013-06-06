% % % Error Results NER
% clear; close all; clc;
% nip_init();
% 
% % Numero de dipolos a considerar
% Nd = 4000; % 1000, 2000, 4000, 8000
% 
% % Cargar datos(lead field, mesh del cerebro etc...
% load(strcat('data/montreal',num2str(Nd),'_10-10.mat'))
% 
% cfg.L = L;
% cfg.cortex = cortex_mesh;
% cfg.fs = 500; % Frecuencia de muestreo para la simulacion
% cfg.t = 0:1/cfg.fs:0.250; % Vector de tiempo
% cfg.elec = elec;
% 
% % Crea una estructura (model) con los datos cargados y declarados arriba
% model = nip_create_model(cfg);
% clear cfg L cortex_mesh eeg_std head elec
% 
% 
% 
% SNR = [-10 0 10];
% NDIP = [1,2,3,4];
% % SNR = [10];
% % NDIP = [2];
% NITER = 30;
% dir = '/home/jscastanoc/Desktop/Results_NER_1010/';
% 
% methods = {'BMF','BMF-TF','ARD','ARD-TF'};
% start = tic;
% for ndip = NDIP
%     for iter = 1:NITER
%         m= 1;
%         tic
%         for snr = SNR
%             %             methods = {'S-FLEX-TF'}
%             for m_idx = 1:length(methods);
%                 method = methods{m_idx};
%                 file_name = strcat(dir,method,'_',num2str(snr),'_',num2str(ndip),...
%                     '_',num2str(iter));
%                 load(file_name);
%                 mse = nip_error_mse(J,J_est);
%                 sai = nip_error_sai(model.cortex,J,J_est,3.5);
%                 tai = nip_error_tai(y,model.L,J_est);
%                 error_tab{1,m_idx}(iter,m) = sai;
%                 error_tab{2,m_idx}(iter,m) = tai;
%                 error_tab{3,m_idx}(iter,m) = mse;
%             end
%             m = m+1;
%         end
%         iter
%         toc
%     end
%     file_name = strcat('error1010_',num2str(ndip),'Nd');
%     save(file_name, 'error_tab')
% end
% toc(start)

%% Plot error
close all;
ndip =4;
file_name = strcat('error1010_',num2str(ndip),'Nd');
load(file_name)

colors ={'-k','--k',':k','-.k','-k'};
% colors = {'b','r','g','k'};
methods = {'BMF','BMF-TF','ARD','ARD-TF'}

Errors = [1,2,3];

for i = Errors
    if i == 3
        fig(i) = figure('Units','normalized','position',[0.2 0.2 0.22 0.2]);
    else
        fig(i) = figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
    end
end

for i = Errors
    for m_idx = 1:length(methods)
        [row,col] = find(isnan(error_tab{i,m_idx}));
        nanMean = nanmean(error_tab{i,m_idx},1);
        if ~isempty(row) && ~isempty(col)
            for m = row'
                for n = col'
                    try
                        error_tab{i,m_idx}(m,n)=nanMean(n);
                    catch
                        disp('hola')
                    end
                end
            end
        end
        error_mean{i,m_idx} = mean(error_tab{i,m_idx});
        X = [1 2 3];
        Y = error_mean{i,m_idx};
        E = std(error_tab{i,m_idx});
        hold on
        %         errorbar(X,Y,E)
    end
    
    figure(fig(i))
    bar(reshape([error_mean{i,:}],3,4))    
    colormap gray
end

titles = {'SAI','TAI','MSE'}
sf = false;
for i = Errors
    figure(fig(i))
%     title(titles{i})
    if i == 1|| i==2
        ylim([0 1])
    end
    if i == 3
        legend(methods,'Location','EastOutside')
    end
    set(gca,'Xtick',[1 2 3])
    set(gca,'Xticklabels',[-10 0 10])
    xlabel('dB')
    ylabel('a.u.')
    if sf
        savefig(strcat(titles{i},num2str(ndip)),gcf,'eps')
    end
end