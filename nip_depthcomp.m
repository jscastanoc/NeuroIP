function varargout = nip_depthcomp(L,varargin)
% varargout = nip_depthcomp(L,varargin)
% Compensate the depth bias on the reconstruction by normalizing the lead
% field matrix.
%  Input:
%       L       -> Ncx3Nd. Lead field matrix.
%       options -> struct
%				type: 'Lnorm' or 'sLORETA'
%				gamma: (only if Lnorm) strengh of normalization
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
        Winv = zeros(1,Nd);
        index = (1:3:Nd);
        for i = index
            norm_term = sqrt((norm(L(:,i),2)^2+ ...
                norm(L(:,i+1),2)^2 + norm(L(:,i+2),2)^2)^gamma);
            Winv(i:i+2) = 1/norm_term;
        end
        Winv = spdiags(Winv', 0, Nd, Nd);
        Winv = Winv / norm(Winv, 'fro');
        Lcomp = full(L*Winv);
    case 'sLORETA'
        Lsloreta = nip_translf(L);
        Winv = sloreta_invweights(Lsloreta);
        inds = reshape(reshape(1:Nd, Nd/3, 3)', [], 1);
        Winv = Winv(inds, inds);
        Winv = Winv / norm(Winv, 'fro');
        Lcomp = full(L*Winv);      
end
if nargout == 2
    extras.Winv = Winv;
    varargout{2} = extras;
end
varargout{1} = Lcomp;
