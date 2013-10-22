clear
% Reduce Kronecker test
L = randn(5,10);
x = eye(10);

% Approach 1
tic
A = [] ;
nonz_idx = find(x);
parfor i = 1:5
    auxL = kron(L(i,:),L);
    auxL = auxL(:,nonz_idx);
    A = [A; auxL];
end
toc

kronLL = kron(L,L);
y1 = kron(L,L)*x(:);
y2 = A*diag(x);
fprintf('Error:%d\n',norm(y1-y2));

% Approach 2
tic
A = [] ;
i = 1;
parfor i = 1:10
    auxL = kron(L(:,i),L(:,i));
    A = [A auxL];
end
toc

kronLL = kron(L,L);
y1 = kron(L,L)*x(:);
y2 = A*diag(x);
fprintf('Error:%d\n',norm(y1-y2));