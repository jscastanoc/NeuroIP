function Q = BMFUr(yp,L)
% Q = BMFUr(yp,L)
%
% Juan S. Castano C.
% 14 Jun 2013
[Nc, Nd] = size(L);
Nr = size(yp,2)
Q = zeros(Nd,Nr);
for i = 1:Nr
    Q(:,i)  = nip_lcmv(yp(:,i),L);
end