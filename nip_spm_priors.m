function [h,Qp] = nip_spm_priors(y,L,patches,Qe,type)
% function h = nip_spm_priors(y,L,patches,Qe,type)
% Find the weights for a set of cortical "patches" using ARD or GS as
% implemented in the SPM software package
% Input:
%       y -> NcxNt. Matrix containing the data,
%       L -> NcxNd. Lead Field matrix
%       Patches -> NdxNp. Each column of the matrix contains the each of
%               the spatial pattern.
%       Qe -> NcxNc. Covariance matrix at the sensor level.
%       type -> string. Method used to find the hyperparameters:
%                   'ARD' - 'GS'
% Output:
%       h -> Npx1. Weights for each of the patches.
%
% Additional Comments: This function is based on a script written by
%           	Jose David Lopez - ralph82co@gmail.com
%				Gareth Barnes - g.barnes@fil.ion.ucl.ac.uk
%				Vladimir Litvak - litvak.vladimir@gmail.com
%        This function is a wrapper to some SPM utilities!
%
% Juan S. Castaño
% 9. May 2013


% YY = y*y';
[Nd, Np] = size(patches);
% Nc = size(L,1);
% Qp = {};
% LQpL = {};
% 
% eyeNd = repmat(eye(Nd),[1,1,3]);
% 
% Lnew = zeros(Nc,Nc);
Qp = zeros(3*Nd,Np*3);
for i = 1:Np
%     i
    if size(L,2)~= size(patches,1);
        Ltemp = nip_translf(L);
        for k = 1:3
            for j = 1:3
                if k==j
%                     Lnew(:,:) = Ltemp(:,:,j)*diag(patches(:,i))*Ltemp(:,:,j)';
                    Q(:,1,j) = patches(:,i);
                else
                    Q(:,1,j) = 0.0*ones(Nd,1,1);
                end                
            end
%             LQpL{end+1}.q = Lnew;
            Qp(:,(i-1)*3+k) = nip_translf(permute(Q,[2 1 3]))';
%             Qp{end+1}.q = nip_translf(permute(Q,[2 1 3]))';
            Q = zeros(Nd,1,3);
        end
    else
%         Qp{end + 1}.q   = patches(:,i);
%         LQpL{end + 1}.q = L*diag(patches(:,i))*L';
    end
end
% Np = numel(Qp);
switch(type)
    case {'GS'}
        % Greedy search over MSPs
%         Np    = length(Qp);
%         Q = zeros(3*Nd,numel(Qp));
%         for i = 1:size(Q,2)
%             Q(:,i) = Qp{i}.q;
%         end
        %         Q = sparse(Q);
%         Q = patches;
        % Multivariate Bayes (Here is performed the inversion)
        %------------------------------------------------------------------
        MVB   = spm_mvb(y,L,[],Qp,Qe,16);
        
        % Accumulate empirical priors (New set of patches for the second inversion)
        %------------------------------------------------------------------
        % MVB.cp provides the final weights of the hyperparameters
        %         Qcp           = Q*MVB.cp;
        %         QP{end + 1}   = sum(Qcp.*Q,2);
        %         LQP{end + 1}  = (L*Qcp)*Q';
        %         LQPL{end + 1} = LQP{end}*L';
        h = diag(MVB.cp);
    case {'ARD'}
        % ReML - ARD (Here is performed the inversion)
        %------------------------------------------------------------------
        [Cy,h,Ph,F] = spm_sp_reml(YY,[],[Qe LQpL],1);
        h = h(2:end);
        % Spatial priors (QP)
        %------------------------------------------------------------------
        % h provides the final weights of the hyperparameters
        %         Ne    = length(Qe);
        %         Np    = length(Qp);
        %         hp    = h(Ne + (1:Np));
        % 		qp    = sparse(0);
        %         for i = 1:Np
        %             if hp(i) > max(hp)/128;
        %                 qp  = qp + hp(i)*Qp{i}.q*Qp{i}.q';
        %             end
        %         end
        %
        %         % Accumulate empirical priors (New set of patches for the second inversion)
        %         %------------------------------------------------------------------
        %         QP{end + 1}   = diag(qp);
        %         LQP{end + 1}  = L*qp;
        %         LQPL{end + 1} = LQP{end}*L';
end