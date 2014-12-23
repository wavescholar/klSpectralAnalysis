function [B, C, O] = nnma(phi, psi, zeta, A, k, varargin)
% GENBREG(PHI, PSI, ZETA, A, K, VARARGIN) -- minimize breg(A, BC)
% 
% [B, C, O] = GENBREG(PHI, PSI, ZETA, A, K, VARARGIN)
% 
% [B, C, O] = GENBREG(PHI, PSI, ZETA, A, K, MAXITER)
% 
% PHI -- The elementwise function
% PSI -- Derivative of PHI
% ZETA -- Second Derivative of PHI
% A    -- Input matrix
% K    -- Rank of approximation
% MAXITER -- Max. number of iterations before stopping
% 
% A \approx B*C, O gives the vector of obj. values at each iteration.
% 
% For example,
% An underoptimized version of the Frob. norm NNMA is:
% 
% GENBREG(@(x) x .* x, @(x) 2*x, @(x) 1, A, K, maxit);
% 
% KL-norm NNMA is
% GENBREG(@(x) x .* log(x), 1 + log(x), 1 ./ x, A, K, maxit);
% 
% Burg-divergence is
% GENBREG(@(x) -log(x), -1 ./ x, 2 ./ (x .* x), A, K, maxit);
% 
% 
% (c) 2005,2006 Suvrit Sra 
% The University of Texas at Austin
% Released under the Gnu Public License. Please see GPL for terms.
% THIS SOFTWARE IS UNSUPPORTED. USE AT YOUR OWN RISK.
% No warranties, express or implied.
% 

[m, n]=size(A);

if (nargin < 6)
  maxiter = 20;
elseif (nargin == 6)
  maxiter = varargin{1};
else
  maxiter = varargin{1};
  disp('Ignoring other input params');
end

iter = 0;

B = 0.5*rand(m, k);
C = 0.5*rand(k, n);

obj = bregdiv(phi, psi, A, B*C);
old = obj + 1;

O(1) = obj;
conj = 0;
% You can change the tolerance here from 1e-3 to higher or lower.
while (iter < maxiter && abs(old-obj) > 1e-3)
  fprintf('Iter: %d %f %f\n', iter, obj,conj);
  a = A(:,1); c = C(:,1);
  t = B*C;
  t1 = zeta(B*C);
  t2 = t1 .* A;
  nr = t2*C';
  
  t3 = t1 .* t;
  dr = t3*C';
  dr = dr + 1e-6;
  m  = nr ./ dr;
  B = B .* m;
  
  t = B*C;
  t1 = zeta(B*C);
  t2 = t1 .* A;
  nr = B'*t2;
  
  t3 = t1 .* t;
  dr = B'*t3;
  dr = dr + 1e-6;
  C = C .* nr ./ dr;
  d = C(:,1);

  obj  = bregdiv(phi, psi, A, B*C); 
  O(iter+1)=obj;
  iter = iter + 1;
end

