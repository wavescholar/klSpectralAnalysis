function f = wishartedf(c)
%WISHARTEDF Density of eigenvalues of an infinite Wishart matrix.
%     THWISHART(C) returns the theoretical limiting density of the 
%                  eigenvalues of an infinite Wishart matrix with 
%                  C = N / M where N and M are the parameters of 
%                  the Wishart matrix.
%
%     C is a positive real number less than 1. 
%     F is a symbolic variable.
%
%     References:
%     [1] Alan Edelman, Handout 3: Experiments with Classical
%                       Matrix Ensembles,  Fall 2004,
%                       Course Notes 18.338.
%     [2] Alan Edelman, Random Matrix Eigenvalues.
%     [3] J. Wishart,   The generalized product moment distribution in
%                       samples from a normal multivariate population,
%                       Biometrika, vol. 20-A, pages 32-52, 1928.
%
%     Alan Edelman and Raj Rao, Sept. 2004.
%     $Revision: 1.0 $  $Date: 2004/09/11  20:49:20 $

syms f x

a1 = (1-sqrt(c))^2;
a2 = (1+sqrt(c))^2;

f = sqrt((x-a1)*(a2-x))/(2*pi*x*c);
