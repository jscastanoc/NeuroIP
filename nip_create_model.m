function model = nip_create_model(cfg)
% model = nip_create_model(cfg) 
% returns a copy of the input structure.
% If the input structure constains L (Leadfield matrix) and t (time vector),
% the returned structure will also contain number of dipoles in L (3Nd), number of 
% channels, and number of time instants.

model = cfg;

if isfield(model,'L')
    model.Nd = size(model.L,2);
    model.Nc = size(model.L,1);
end
if isfield(model,'t')
    model.Nt = length(model.t);
end


