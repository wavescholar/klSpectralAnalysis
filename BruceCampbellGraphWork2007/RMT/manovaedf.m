function f = manovaedf(c1,c2);
%MANOVAEDF Theoretical density of eigenvalues of an infinite MANOVA matrix.
%     MANOVAEDF(C1,C2) returns the limiting theoretical density of the  
%                      eigenvalues of an infinite MANOVA matrix,
%                      where C1 = N / M1 and C2 = N / M2 and N1, M1, and 
%                      M2 are parameters of the MANOVA matrix. 
%
%     C1 and C2 are positive real numbers less than 1.
%     F is a symbolic variable.
%
%     References:
%     [1] Alan Edelman,   Handout 3: Experiments with Classical
%                         Matrix Ensembles, Fall 2004,
%                         Course Notes 18.338.
%     [2] Alan Edelman,   Random Matrix Eigenvalues.
%     [3] R. J. Muirhead, Aspects of Multivariate Statistical Theory,
%                         John Wiley & Sons, New York, 1982.
%
%
%     Alan Edelman and Raj Rao, Sept. 2004.
%     $Revision: 1.0 $  $Date: 2004/09/10  23:55:18 $

syms x f

b0 = c1*x-c2*x-c1+2;
b1 = -2*c2*x^2+2*x-3*c1*x+c1+c2*x-1+2*c1*x^2;
b2 = c1*x-2*c1*x^2+c2*x^2-x^3*c2+x^3*c1;

f = sqrt(4*b2*b0-b1^2)/(2*pi*b2);
