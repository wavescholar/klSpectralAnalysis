% SIGMA_CALC_WEIGHTS
% Ross M. Richardson
% Fan Chung Graham's Research Group, UCSD
% Begun 7/13/05

function [Weights,Laplacian,sigma] = ...
    sigma_calc_weights(G, epsilon)

% This file is part of a private research program to compute
% the sigma function and associated weight vector of a graph G
% (see Spectral Graph Theory, by Dr. Fan Chung Graham, pg. 104).
%
% Input: G - nxn 0/1 (sparse) adjacency matrix of a graph G.
% Only the upper triangular portion is required.  Note that
% the diagonal entries should be 0.
%
% Output: W - nxn weight matrix , w_ij >= 0
%         sigma - sigma value for G

[A,b,c,K,slacks] = sigma_make_instance(G);
[x,y,info] = sedumi(A,b,c,K);

sigma = -c'*x;

[m,n] = size(x);

% W is the PSD matrix which optimizes our SDP
W=mat(x(slacks+1:m,1));

% F is the positive square root of the diagonal of W
F = sqrt(diag(diag(W)));

temp = inv(F) * W * inv(F);
lambda = -nthroot(temp(1,1), 3);

% Now we get a matrix which is a_ij >= 0 only if i ~ j
A = (temp - diag(diag(temp))) / (lambda * lambda);

% Normalize A by the largest eigenvalue
B = A / max(eigs(A));

% Now we find a non-negative fixed vector of B and use it to
% construct the weight matrix.
[m,n] = size(B);
[V,D] = eigs(B);
[mm,nn] = size(V);

v = zeros(n,1);
for i = 1:nn
  if(D(i,i) >= 1.0)
    % Fixed vector.  Now test to see if all positive.
    v = V(:,i);
    if (sum(v + abs(v)) == 0 | sum(v - abs(v)) == 0)
      break % Found our vector
    end
  end
end

v = abs(v);

for i = 1:n
  for j = 1:n
    if (i == j)
      Weights(i,j) = 0;
    else
      Weights(i,j) = B(i,j) * v(i) * v(j);
    end
  end
end

for i = 1:n
  for j = 1:n
    if (i == j)
      Laplacian(i,j) = 1;
    else
      temp = sqrt(sum(Weights(i,:)) * sum(Weights(j,:)));
      if (abs(temp) > epsilon)
        Laplacian(i,j) = - Weights(i,j) / temp;
      else
        Laplacian(i,j) = 0;
      end
    end
  end
end

% Uncomment to see "Laplacian" I - B
%  eye(n) - B
