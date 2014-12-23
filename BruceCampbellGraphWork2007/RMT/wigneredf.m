function f = wigneredf;
%WIGNEREDF Density of eigenvalues of an infinite Wigner matrix.
%     WIGNEREDF returns the theoretical density of the eigenvalues of 
%               an infinite Wigner matrix
%
%     F is a symbolic variable.
%
%     References:
%     [1] Alan Edelman, Handout 3: Experiments with Classical
%                       Matrix Ensembles, Fall 2004,
%                       Course Notes 18.338.
%     [2] Alan Edelman, Random Matrix Eigenvalues.
%     [3] E. P. Wigner, Characteristic vectors of bordered matrices with
%                       infinite dimensions, Annals of Mathematics,
%                       vol. 62, pages 548-564, 1955.
%
%     Alan Edelman and Raj Rao, Sept. 2004.
%     $Revision: 1.0 $  $Date: 2004/09/11  20:47:15 $

syms f x

f = sqrt(4-x^2)/(2*pi);
