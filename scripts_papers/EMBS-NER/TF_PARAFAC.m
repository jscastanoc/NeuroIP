function Factors = TF_PARAFAC(y)
% Factors = TF_PARAFAC(y)
%
% Juan S. Castano
% 25 May 2013


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time Frequency decomposition %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 5; % Time shift
M = 250; % Frequency Res
c = dgtreal(y','gauss',a,M);
c = permute(c,[3 1 2]);

% %%%%%%%%%%%
% % Show TF %
% %%%%%%%%%%%
% ch = 30; %% Show the TF of this channel
% figure('Units','Normalized','Position',[0.3 0.1 0.3 0.3])
% subplot(2,1,1)
% plot(model.t,model.y(ch,:))
% subplot(2,1,2)
% % c_sgram=sgram(model.y(ch,:),model.fs,'lin','wlen',30);
% cc(:,:) = c(ch,:,:);
% plotdgtreal(cc,a,M,'linsq');



% Parafac Decomposition
Opt(1) = 1e-6; Opt(2) = 1; Opt(3) = 0; Opt(4) = 0; Opt(5) = 10; Opt(6) = 2500;
const = [2,2,2];

Nfac =2; % Number of factors to decompose
concord = 100;
while concord > 95
    [Factors,it,err,concord] = parafac(abs(c),Nfac,Opt,const);
    if concord <95 && Nfac == 2
        Nfac = 1;
        [Factors,it,err,concord] = parafac(abs(c),Nfac,Opt,const);
        break
    end
    Nfac = Nfac+1;
end
% %%%%%%%%%%%%%%%%%%%%%%%
% % Show PARAFAC DECOMP %
% %%%%%%%%%%%%%%%%%%%%%%%
% figure('Units','Normalized','Position',[0.3 0.1 0.3 0.3])
% title('Factors')
% subplot(1,3,1)
% plot(Factors{1})
% subplot(1,3,2)
% plot(Factors{2})
% subplot(1,3,3)
% plot(Factors{3})

end