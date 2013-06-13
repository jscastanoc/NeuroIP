function promCor = nip_error_mse(Js,Jr)
% promCor = nip_error_mse(Jsim,Jrec)
% 
% Juan S. Castano C
% 22 May 2013

% Js = Jsim(find(Mr),:);
% Jr = Jrec(find(Mr),:);

% recV = Jr*Js';
% recV = norm(sum((Js-Jr).^2,2));

% realV = Js*Js';
% realV = norm(sum(Js.^2,2));

% promCor = sum(diag(recV).^2)/sum(diag(realV).^2);
promCor = mean(sum((Jr-Js).^2,2));
% promCor = recV/realV;
