clear all
close all

A = rand(3,3);
A = A*A';
A = A - diag(diag(A));
D = sum(A,2);
DD = diag(D);
K = DD - A;

sD = D.^0.5;
iD = sD.^-1;

K - diag(D) + A
diag(iD)*diag(D)*diag(iD)

L = diag(iD) * A * diag(iD);

%[u1,s1,v1] = svd(L)
[u1,s1] = eig(L)

%[u2,s2,v2] = svd(K)
[u2,s2] = eig(K)

diag(iD)*K*diag(iD) - eye(3) + u1*s1*v1'


u3 = diag(sqrtD)*u1;
u3*(eye(3)-s1)*u3'



