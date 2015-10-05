function Ar = markovSuppress(A,trs)
  
  
  A = sparse(A);
  D = sum(A,1)';                 % Normalize column sum to one.
  R = A * sparse(diag(D.^-1));  % Markov

  Tol = trs*sum(R+R',2) * ones(1,size(R,2));
  RR  = (R + R') > 2*Tol;
  % find locations where p(1->2) is 0 but not p(2->1)
  RR(find(abs(RR-RR'))) = 1;
  
  Ar = A .* RR ;
