function [U,S] = powIter(L,U)
% power iteration to update eigenvectors.
% orthogonalize with fixed gram procedure.
% the number of power iterations = 2*size(U,2)
% the eigen values are obtaineed by rayleigh coefficient
% measure. can the rayleigh coefficients be ever negative?
%
% to do: stopping criterion for power iterations
% calls: gramFixed (to orthogonalize power iterated basis)
%
  tt = size(U,2);
  for k = 1:tt %2*tt
    U = L*U;
  end
  
  [id,U] = gramFixed(U,0);
  S = diag(U'*L*U);
  [S,id] = sort(-S);
  S = -S;
  U = U(:,id);
