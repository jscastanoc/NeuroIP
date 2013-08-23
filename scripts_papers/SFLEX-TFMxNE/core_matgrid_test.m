function [J_rec time] = core_matgrid_test(model,method)

iter_basis = [1.5];

switch method
    case 'LOR'
        %                     [Laplacian] = nip_neighbor_mat(model.cortex);
        %                     Q = inv(Laplacian*Laplacian');
        Q = speye(model.Nd);
        [J_rec,~] = nip_loreta(model.y,model.L,Q);
        
    case 'S-FLEX'
        [S, out] = sflex_cortical_dal(model.y, nip_translf(model.L),model.cortex,...
            struct('eps',0.1,'sigma',iter_basis));
        J_rec = nip_trans_solution(S);
        
    case 'TF-MxNE'
        options.iter = 50;
        options.spatial_reg = 0.9;
        options.temp_reg = 0.15;
        options.tol = 1e-2;
        options.a = 10;
        options.M = 200;
        [J_rec,~] = nip_tfmxne_port(model.y,model.L,options);
        
    case 'S+T'
        options.iter = 50;
        options.spatial_reg = 0.9;
        options.temp_reg = 0.15;
        options.tol = 1e-2;
        options.a = 10;
        options.M = 200;
        [J_rec,~] = nip_sflex_tfmxne(model.y,model.L,basis{1},options);
    otherwise
        error(strcat('Nah! ',method,' is not available'))
end
time = toc;
end