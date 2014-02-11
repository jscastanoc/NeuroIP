function data_m = nip_energy(data)
% function data_m = nip_energy(data)
% Auxiliary function that returns the magnitude of each dipole vector.
% Input:
%       data -> 3NdxNt. Vector containing the activity of Nd dipoles at a
%       given time sample. The format or ordering of the data is.
%               (x1,y1,z1,....xNd,yNd,zNd)'. where xn yn and zn is the activity
%               of the n-th dipole in each coordinate.
% Output:
%       data_m -> Ndx1. Magnitude of the activity in each dipole.
%
% Juan S. Castano C.
% jscastanoc@gmail.com
% 2 Sep 2013.
[Nd Nt] = size(data);
data_m = zeros(Nd/3,1);
for j = 1:Nt
for i = 1:Nd/3
    data_m(i,j) = sqrt(sum(data((i-1)*3+1:(i-1)*3+3,j).^2));
end
end
