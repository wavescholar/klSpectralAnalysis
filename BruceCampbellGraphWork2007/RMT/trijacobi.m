function t = trijacobi(n,a,b)
%TRIJACOBI Returns the tridiagonal matrix of polynomial recurrence coefficients
%     TRIJACOBI(N,A,B) returns a tridiagonal matrix whose entries correspond  
%                      to the coefficients describing the recurrence between
%                      the Jacobi Polynomials. 
%
%     The eigenvalues of this matrix correspond to the roots of the N-th        
%     Jacobi polynomial with parameters A and B.       
%                       
%     N, A and B are integers.
%     
%
%     References:
%     [1] Alan Edelman,   Handout 6: Tridiagonal Matrices, Orthogonal 
%                         Polynomials ans the Classical Random Matrix 
%                         Ensembles, Fall 2004,
%                         Course Notes 18.338.
%     [2] Alan Edelman,   Random Matrix Eigenvalues.
%     [3] Gabor Szego,    Orthogonal Polynomials, American Mathematical
%                         Society, Providence, 1975. 4th Edition.
%     [4] Eric Weisstein, Jacobi Polynomial." From MathWorld--
%                         A Wolfram Web Resource. 
%                         http://mathworld.wolfram.com/JacobiPolynomial.htm
%
%     Brian Sutton, Sept. 2004.
%     $Revision: 1.0 $  $Date: 2004/09/23  11:35:18 $

i = (1:n-1)';

d1 = 2*sqrt(i.*(a+i).*(b+i).*(a+b+i)./(-1+a+b+2*i)./(1+a+b+2*i))./(a+ ...
                          b+2*i);
d0 = -(a^2-b^2)./(a+b+2*(0:n-1)')./(2+a+b+2*(0:n-1)');

t = spdiags([[d1;0] d0 [0;d1]],[-1 0 1],n,n);
