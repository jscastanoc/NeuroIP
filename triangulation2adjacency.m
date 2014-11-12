function A = triangulation2adjacency(face,vertex)

% triangulation2adjacency - compute the adjacency matrix
%   of a given triangulation.
%
%   A = triangulation2adjacency(face);
% or for getting a weighted graph
%   A = triangulation2adjacency(face,vertex);
%
%   Copyright (c) 2005 Gabriel Peyr�
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

% [tmp,face] = check_face_vertex([],face);
f = double(face)';

A = sparse([f(:,1); f(:,1); f(:,2); f(:,2); f(:,3); f(:,3)], ...
           [f(:,2); f(:,3); f(:,1); f(:,3); f(:,1); f(:,2)], ...
           1.0);
% avoid double links
A = double(A>0);

return; 


nvert = max(max(face));
nface = size(face,1);
A = spalloc(nvert,nvert,3*nface);

for i=1:nface
    for k=1:3
        kk = mod(k,3)+1;
        if nargin<2
            A(face(i,k),face(i,kk)) = 1;
        else
            v = vertex(:,face(i,k))-vertex(:,face(i,kk));
            A(face(i,k),face(i,kk)) = sqrt( sum(v.^2) );    % euclidean distance
        end
    end
end 
% make sure that all edges are symmetric
A = max(A,A');
