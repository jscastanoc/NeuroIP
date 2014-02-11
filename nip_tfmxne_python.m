function [J_rec extras] = nip_tfmxne_python(y,L,varargin)
%NIP_TFMXNE_PYTHON Summary of this function goes here
%   Detailed explanation goes here
% pycode = fileread('test.py');
% py_export('y','L','alpha_space', 'alpha_time');
% py('debugon')
% py('eval',pycode);
% X = py('get','X');
% active_set = py('get','active_set');

Lor = L;
yor = y;

p = inputParser;

def_sreg = 80;
def_treg= 1;
def_maxiter = 300;
def_tol = 1e-6;
def_resnorm = 0.3;
def_lipschitz = [];
def_optimres = false;
def_Winv = [];
def_wsize = double(64.0);
def_tstep = double(8.0);

addParamValue(p,'tstep',def_tstep);
addParamValue(p,'wsize',def_wsize);
addParamValue(p,'sreg',def_sreg);
addParamValue(p,'treg',def_treg);
addParamValue(p,'maxiter',def_maxiter);
addParamValue(p,'tol',def_tol);
addParamValue(p,'resnorm',def_resnorm);
addParamValue(p,'lipschitz',def_lipschitz);
addParamValue(p,'optimres',def_optimres);
addParamValue(p,'Winv',def_Winv)

parse(p,varargin{:})
options = p.Results;


tfmxne = pyimport('tfmxne');
numpy = pyimport('numpy');

[Nc Nd] = size(L);
Nt = size(y,2);
y = numpy.asarray(y);
L = numpy.asarray(L);

sreg = options.sreg;
treg = options.treg;
options.tol
rev_line = '';
for i = 1:20% Adjust the regularization parameters 20 times max.
    clck = tic;
    out = tfmxne.tfmxne(y,L,sreg,treg,...
        options.wsize,options.tstep,options.maxiter,options.tol);
    
    out =  unpy(out);
    
    X = out{1};
    active_set = out{2};
    J_rec = zeros(Nd,Nt);
    J_rec(find(active_set),:) = X;
    extras.active_set = active_set;
    
      
    
    J_rec = J_rec*(norm(yor,'fro')/norm(Lor*J_rec,'fro'));
    resnorm = norm(yor-Lor*J_rec, 'fro')/norm(yor, 'fro');
    
    msg = sprintf('Iteration # %d, Time per iter: %f \nDipoles~=0: %d \nRegPar S:%d T:%d \nResidual Norm = %8.5e\n'...
        ,i,toc(clck),sum(full(active_set))/3, sreg,treg,resnorm);
    fprintf([rev_line, msg]);    
    rev_line = repmat(sprintf('\b'),1,length(msg));
    
    if (resnorm < options.resnorm || ~options.optimres)
        break;
    else
        sreg = 0.7*sreg;
        treg = 0.7*treg;
    end
end

if ~isempty(options.Winv)
    J_est = nip_translf(J_rec');
    siJ = size(J_est);
    J_rec = permute(reshape(full(reshape(permute(J_est, [1 3 2]), siJ(1), [])*options.Winv), siJ(1), siJ(3), siJ(2)), [1 3 2]);
    J_rec = nip_translf(J_rec)';
end

end

