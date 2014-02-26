function [J_rec, extras] = nip_sflex(y, L, basis, varargin)
%  [J_est, extras] = nip_sflex(y, L, basis, reg_par)
% Implements "Large-scale EEG/MEG source localization with spatial
% flexibility." by Haufe et al 2011
%
% Input:
%       y -> NcxNt. Matrix containing the data,
%       L -> NcxNd. Lead Field matrix
%       basis-> NdxNs. Matrix containing the spatial basis functions.
%       reg_par-> Scalar. Regularization parameter (1e-2 by default)
% Output:
%       J_rec -> NdxNt. Reconstructed activity (solution)
%       extras.regpar -> Scalar. Optimum regularization parameter
%
% Additional comments: Uses the DAL optimization toolbox.
%
% Juan S. Castano C.
% 13 June 2013
NDUM = 3;
yor = y;
Lor = L;

p = inputParser;
def_maxiter = 50;
def_resnorm = 0.5;
def_regpar = 100;
def_optimres = false;
def_Winv = [];
addParamValue(p,'maxiter',def_maxiter);
addParamValue(p,'resnorm',def_resnorm);
addParamValue(p,'regpar',def_regpar);
addParamValue(p,'optimres',def_optimres);
addParamValue(p,'Winv',def_Winv);

parse(p,varargin{:})
options = p.Results;

[Nc Nt] = size(y);
Nd = size(L,2);

if nargin <=3
    reg_par = 1e-2;
end

L = nip_translf(L);
for i = 1:3
    L(:,:,i) = L(:,:,i)*basis; % J simulado FINAL
end
nbasis = size(basis,2);

L = nip_translf(L);

% A = sparse(kron(speye(Nt), L));

% A = nip_translf(A);
% A = permute(A,[1 3 2]);



% [xx0,~] = nip_loreta(y,L,diag(nip_lcmv(y,L)));
% xx0f =reshape(xx0',[3*Nt,nbasis]);

xx0f = zeros(3*Nt,nbasis);
% [xx,status]=dalsqgl(zeros(3,nbasis*Nt), A, y(:), reg_par);

XX = {@xforth, @xback, Nc*Nt, nbasis*NDUM*Nt};

opt.solver = 'qn';
opt.maxiter = 1000;
opt.tol = 1e-8;
%opt.stopcond = 'fval';
iter = 1;
xx = xx0f;
while true
    [xx,status]=dalsqgl(xx, XX, y(:), options.regpar, opt);
    % xx = xx(:);
    % xx = reshape(xx,[nbasis*3,Nt]);    
    index = 1:3:Nd;
    clear xxf;
    for i = 0:2
        xxf(index+i,:) = xx(1+i:3:end,:)';
        xxf(index+i,:) = basis*xxf(index+i,:); % J simulado FINAL
    end
    J_rec = xxf*(norm(yor,'fro')/norm(Lor*xxf,'fro'));    
    resnorm = norm(yor-Lor*J_rec, 'fro')/norm(yor, 'fro');
    fprintf('GOF = %8.5e\n', resnorm);
    
    if iter > options.maxiter || resnorm < options.resnorm || ~options.optimres
        break;
    else
        options.regpar = 0.75*options.regpar;
    end
end

if ~isempty(options.Winv)
    J_est = nip_translf(J_rec');
    siJ = size(J_est);
    J_rec = permute(reshape(full(reshape(permute(J_est, [1 3 2]), siJ(1), [])*options.Winv), siJ(1), siJ(3), siJ(2)), [1 3 2]);
    J_rec = nip_translf(J_rec)';
end
extras =status;

    function xfo = xforth(x)
        Q = size(x, 2);
        [in1 indum] = find(x);
        in2 = unique(ceil(in1./(Nt)));
        l = length(in2)/NDUM;
        if l == 0
            xfo = zeros(Nc*Nt, Q);
        else
            xfo = reshape((L(:, in2)*reshape(permute(reshape(full(x(in1, :)), NDUM, Nt, l, Q), [1 3 2 4]), NDUM*l, Nt*Q)), Nc*Nt, Q);
        end
    end

    function xba = xback(x)
        xba = reshape(permute(reshape(L'*reshape(x, Nc, Nt), NDUM, nbasis, Nt), [1 3 2]), NDUM*Nt*nbasis, []);
    end
end