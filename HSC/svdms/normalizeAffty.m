function [A,D,L,M,H] = normalizeAffty(A,medFlag) % H added to output by virginia, 09-24
%function [A,D,L,M] = normalizeAffty(A)
% returns sparse versions of matrices
%  A : sparse affinity matrix  (normalizes by the median)
%  D : sum(A,1)' as a column vector
%  L : D^-0.5*A*D^-0.5 (normalized affty)
%  M : A*D^-1 (markov)
%  H : Kirchoff matrix D - A
% 
  if ~issparse(A)
    A = sparse(A);
  end
  if nargin > 1 
    if medFlag
    A = A/median(full(sum(A,1)));
    %A = A/sum(A(:));
    end
  end
  
  D = full(sum(A,1))';                 % Normalize column sum to one.
  M = A * spdiags(D.^-1, 0, length(D), length(D));  % Markov
  %sqrtD    = spdiags(D.^0.5, 0 , length(D),length(D));
  sqrtDinv = spdiags(D .^ -0.5, 0, length(D),length(D) );
  L = sqrtDinv * A * sqrtDinv;
  H = spdiags(D,0,length(D),length(D)) - A; % brought back by virginia, 09-24

return;


