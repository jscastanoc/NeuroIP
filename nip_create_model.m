function model = nip_create_model(cfg)
% model = nip_create_model(cfg) 

model = cfg;

if isfield(model,'L')
    model.Nd = size(model.L,2);
    model.Nc = size(model.L,1);
end
if isfield(model,'t')
    model.Nt = length(model.t);
end


