function er = nip_all_errors(y,L,J_rec,Jclean,cortex,actidx)
% er = nip_all_errors(y,L,J_rec,Jclean,cortex,actidx)
% Function to compute Earth mover's distance, maximum correlation, geodesic distance,
% weighted geodesic distance, explained variance and spatial accuracy inde.
% Input:
% 		y -> NcxNt. Matrix with EEG measurements.
% 		L -> 3NdxNt. Lead field matrix.
% 		J_rec -> 3NdxNt. Reconstructed activity.
% 		Jclean -> 3NdxNt. Ground truth of the brain activity.
% 		cortex -> Struct. Structure containing the mesh describing the cortex of the brain
% 						 With the fields vc and tri or vertices and faces.
% 		actidx -> Indices of the active dipoles (between 1 and Nd).
%
% Output:
% 		er -> 1x9 vector with the following errors:
%             1) EMD
%             2) Max. Correlation
%             3) Average distance between local maxima X
%             4) 3. weighted with the achieved correlation X
%             5) Explained variance
%             6) Spatial accuracy
%             7) Specificity
%             8) Sensitivity
%             9) Root mean squared error
%
% Juan S. Castano
% jscastanoc@gmail.com

Nd = size(cortex.vc,1);
tic
distmat = nip_fuzzy_sources(cortex,[],struct('dataset','montreal','save',true,'calc','dist'));
toc
er = [];


actsources = length(actidx);

sig1 = sqrt(sum(J_rec.^2,2));
sig1 = sig1/norm(sig1);
sig2 = sqrt(sum(Jclean.^2,2));
sig2 = sig2/norm(sig2);
er(1) = nip_emd(sig1,sig2,distmat);
Jrecbckp = J_rec;
J_rec = J_rec/norm(J_rec,2);
Jclean = Jclean/norm(Jclean,2);

CR = Jrecbckp;
for i = 1:actsources
    CR = [Jclean((actidx(i)-1)*3+1:(actidx(i)-1)*3+3,:); CR];
end
CR = corr(CR');
CR = CR - eye(length(diag(CR)));
[val, idxs] = max(abs(CR(1:actsources*3,actsources*3+1:end)),[],2);
er(2) = nanmean(val);

distmax = 0;
weighted_d = 0;
for i = 1:actsources*3
    hhdist = distmat(ceil(idxs(i)/3),actidx(ceil(i/3)));
    if val(i) == 0
        aux = NaN;
    else
        aux = 1/val(i);
    end
    weighted_d(i) = hhdist*aux;
    distmax(i) = hhdist;
end

er(3) = nanmean(distmax);
er(4) = nanmean(weighted_d);
er(5) = nip_error_tai(y,L,Jrecbckp);
[out,~,~] = nip_error_sai(cortex, Jclean,J_rec,6);
er(6) = out.sai;
er(7) = out.spec;
er(8) = out.sens;
out = nip_calc_error(y,L,Jclean,Jrecbckp);
er(9) = out.rmse;