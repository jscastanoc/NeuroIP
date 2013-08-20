function distance = nip_emd(Sol1, Sol2,cortex)

if ~isfield(cortex, 'vc') && ~isfield(cortex,'tri')
    cortex.vc = cortex.vertices;
    cortex.vertices = [];
    cortex.tri = cortex.faces;
    cortex.faces = [];
end


Nd = length(Sol1);

data_m = zeros(Nd/3,1);
for i = 1:Nd/3
    data_m(i) = sqrt(sum(Sol1((i-1)*3+1:(i-1)*3+3).^2));
end
sig1 = data_m;

data_m = zeros(Nd/3,1);
for i = 1:Nd/3
    data_m(i) = sqrt(sum(Sol2((i-1)*3+1:(i-1)*3+3).^2));
end
sig2 = data_m;

fprintf('Computing Earth s moving distance... ')
tic
aff = graphrbf(cortex);
% aff = 1.15.^(graphrbf(cortex));
% aff = aff - diag(diag(aff));
distance= emd_hat_mex(sig1,sig2,aff);
% distance= emd_mex(sig1',sig2',aff);
fprintf('done! \nElapsed time: %.2d secs \n', toc)
end