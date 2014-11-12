function W = compute_mesh_weight(vertex,face,type,options)

% compute_mesh_weight - compute a weight matrix
%
%   W = compute_mesh_weight(vertex,face,type,options);
%
%   W is sparse weight matrix and W(i,j)=0 is vertex i and vertex j are not
%   connected in the mesh.
%
%   type is either 
%       'combinatorial': W(i,j)=1 is vertex i is conntected to vertex j.
%       'distance': W(i,j) = 1/d_ij^2 where d_ij is distance between vertex
%           i and j.
%       'conformal': W(i,j) = cot(alpha_ij)+cot(beta_ij) where alpha_ij and
%           beta_ij are the adjacent angle to edge (i,j)
%
%   If options.normalize=1, the the rows of W are normalize to sum to 1.
%
%   Copyright (c) 2007 Gabriel Peyre

% Copyright (c) 2009, Gabriel Peyre
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

options.null = 0;
[vertex,face] = check_face_vertex(vertex,face);

nface = size(face,1);
n = max(max(face));

verb = getoptions(options, 'verb', n>5000);

if nargin<3
    type = 'conformal';
end

switch lower(type)
    case 'combinatorial'
        W = triangulation2adjacency(face);
    case 'distance'
        W = my_euclidean_distance(triangulation2adjacency(face),vertex);
        W(W>0) = 1./W(W>0);
        W = (W+W')/2; 
    case 'conformal'
        % conformal laplacian
        W = sparse(n,n);
        ring = compute_vertex_face_ring(face);
        for i = 1:n
            if verb
                progressbar(i,n);
            end
            for b = ring{i}
                % b is a face adjacent to a
                bf = face(:,b);
                % compute complementary vertices
                if bf(1)==i
                    v = bf(2:3);
                elseif bf(2)==i
                    v = bf([1 3]);
                elseif bf(3)==i
                    v = bf(1:2);
                else
                    error('Problem in face ring.');
                end
                j = v(1); k = v(2);
                vi = vertex(:,i);
                vj = vertex(:,j);
                vk = vertex(:,k);
                % angles
                alpha = myangle(vk-vi,vk-vj);
                beta = myangle(vj-vi,vj-vk);
                % add weight
                W(i,j) = W(i,j) + cot( alpha );
                W(i,k) = W(i,k) + cot( beta );
            end
        end
    otherwise
        error('Unknown type.')
end

if isfield(options, 'normalize') && options.normalize==1
    W = diag(sum(W,2).^(-1)) * W;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta = myangle(u,v);

du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
du = max(du,eps); dv = max(dv,eps);
beta = acos( sum(u.*v) / (du*dv) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = my_euclidean_distance(A,vertex)

if size(vertex,1)<size(vertex,2)
    vertex = vertex';
end

[i,j,s] = find(sparse(A));
d = sum( (vertex(i,:) - vertex(j,:)).^2, 2);
W = sparse(i,j,d);  
end

function [vertex,face] = check_face_vertex(vertex,face, options)

% check_face_vertex - check that vertices and faces have the correct size
%
%   [vertex,face] = check_face_vertex(vertex,face);
%
%   Copyright (c) 2007 Gabriel Peyre

vertex = check_size(vertex,2,4);
face = check_size(face,3,4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = check_size(a,vmin,vmax)
if isempty(a)
    return;
end
if size(a,1)>size(a,2)
    a = a';
end
if size(a,1)<3 && size(a,2)==3
    a = a';
end
if size(a,1)<=3 && size(a,2)>=3 && sum(abs(a(:,3)))==0
    % for flat triangles
    a = a';
end
if size(a,1)<vmin ||  size(a,1)>vmax
    error('face or vertex is not of correct size');
end
end