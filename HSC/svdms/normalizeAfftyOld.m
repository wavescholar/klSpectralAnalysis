function [A,D,L,M] = normalizeAfftyOld(A)
%function [A,D,L,M] = normalizeAffty(A)
% returns sparse versions of matrices
%  A : sparse affinity matrix  (normalizes by the median)
%  D : sum(A,1)' as a column vector
%  L : D^-0.5*A*D^-0.5 (normalized affty)
%  M : A*D^-1 (markov)

  if ~issparse(A)
    A = sparse(A);
  end
  A = A/median(full(sum(A,1)));
  %A = A/sum(A(:));

  D = sum(A,1)';                 % Normalize column sum to one.
  M = A * spdiags(D.^-1, 0, length(D), length(D));  % Markov
  %sqrtD    = spdiags(D.^0.5, 0 , length(D),length(D));
  sqrtDinv = spdiags(D .^ -0.5, 0, length(D),length(D) );
  L = sqrtDinv * A * sqrtDinv;
end
