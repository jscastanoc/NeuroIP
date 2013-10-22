function varargout = nip_depthcomp(L,varargin)
% varargout = nip_depthcomp(L,varargin)
% Compensate the depth bias on the reconstruction by normalizing the lead
% field matrix. 
%  Input:
%       L       -> Ncx3Nd. Lead field matrix.
%       gamma   -> Scalar. Normalization parameter. 0 means no
%           normalization. 1 means maximum normalization.
%  Output:
%       Lcomp   -> Ncx3Nd. Normalized lead field matrix.
%       extra
% Juan S. Castano C.
% jscastanoc@gmail.com
% 14 Aug 2013

[Nc Nd] = size(L);

if isempty(varargin);
    varargin{1} = [];
end

if ~isfield(varargin{1},'type')
    vararging{1}.type = 'Lnorm';
end
if ~isfield(varargin{1},'gamma') && strcmp(varargin{1}.type,'Lnorm')
    varargin{1}.gamma = 0.6;    
end


switch varargin{1}.type
    case 'Lnorm'
        gamma = varargin{1}.gamma;
        index = (1:3:Nd);
        for i = index
            norm_term = sqrt((norm(L(:,i),2)^2+ ...
                norm(L(:,i+1),2)^2 + norm(L(:,i+2),2)^2)^gamma);
            Lcomp(:,i:i+2) = L(:,i:i+2)/norm_term;
        end
    case 'sLORETA'
        Lsloreta = nip_translf(L);
        Winv = full(sloreta_invweights(Lsloreta));
        Winv = cat(3,Winv(1:end/3,1:end/3),Winv(end/3 +1:2*end/3,end/3 +1:2*end/3),Winv(2*end/3 +1 :end,2*end/3 +1:end));
        for i = 1:3
            Lsloreta(:,:,i) = Lsloreta(:,:,i)*Winv(:,:,i);
        end
        Lcomp = nip_translf(Lsloreta);        
        if nargout == 2
            extras.Winv = Winv;
            varargout{2} = extras;
        end
end
varargout{1} = Lcomp;