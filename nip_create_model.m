function model = nip_create_model(cfg)
% model = nip_create_model(cfg, type) 

model = cfg;
model.Nd = size(model.L,2);
model.Nc = size(model.L,1);
model.Nt = length(model.t);


