function obj = bregdiv(phi, psi, X, Y)
% 

obj = 0;

[m, n] = size(X);

for i=1:m
  for j=1:n
    x = X(i, j);
    y = Y(i, j);
    obj = obj + (phi(x) - phi(y) - (psi(y)*(x-y)));
  end
end


