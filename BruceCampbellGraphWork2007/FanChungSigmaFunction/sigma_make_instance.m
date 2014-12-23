% SIGMA_MAKE_INSTANCE
% Ross M. Richardson
% Fan Chung Graham's Research Group, UCSD
% Begun 7/12/05

function [A,b,c,K,slacks] = ...
    sigma_make_instance(G)

% This file is part of a private research program to compute
% the sigma function and associated weight vector of a graph G
% (see Spectral Graph Theory, by Dr. Fan Chung Graham, pg. 104).
%
% Input: G - nxn 0/1 (sparse) adjacency matrix of a graph G.
% Note that for space considerations, we use only the upper
% triangular portion of G, e.g. G(1,2) is used, G(2,1) is not.
%
% Output: A SeDuMi instance. slacks = edges + n
%
% Technical Details: The instance sets up the following
% optimization problem.
% min sum -x_ij
% trace X = 1
% x_ij <= 0 if i is adjacent to j or i=j
% X is psd
% See the code for implementation and commentary.

[m,n] = size(G);
if(m ~= n)
  error('Not an nxn matrix.')
end

% Count the number of edges
edges = 0;
for i = 1:n
  for j = i:n
    if(G(i,j) == 1)
      edges = edges + 1;
    end
  end
end

% This cajoles our min problem to a max problem
% c*x will be NEGATIVE the actual value as a result.
% We have edges slack variables and n*n decision variables
% (the psd matrix).
c = -sparse([zeros(edges+n,1); ones(n*n,1)]);

b = sparse([1; zeros(n*(n+1)/2,1)]);

% For A, the first line gives the trace=1 condition
% Then we enforce the x_ij <=0 conditions with the use
% of slack variables.
A = sparse([zeros(1,edges+n), vec(spdiags(ones(n,1),0,n,n))';
            zeros(n*(n+1)/2, (edges + n) + n*n)]);


  
row = 2;
nextslack = 1;
for i = 1:n
  for j = i:n
    if (G(i,j) == 1 | i == j) % Positivity Constraint
      A(row, nextslack) = -1;
      nextslack = nextslack + 1;
    end
    A(row, (edges + n) + (j-1)*n + i) = 1;
    row = row + 1;    
  end
end


% This is test code below--feel free to wipe.
if (1==0)
% Theta instead of sigma
row = 2;
for i = 1:n
  for j = i:n
    if (G(i,j) == 0 & i ~= j)
      A(row, (edges + n) + (j-1)*n + i) = 1;
    end
    row = row + 1;
  end
end

end


% End Theta

    
K.l = [edges+n]; % Number of slack variables
K.s = [n];

slacks = edges + n;
return
