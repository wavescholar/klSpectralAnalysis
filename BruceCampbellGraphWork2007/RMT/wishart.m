function w = wishart(n,m,isreal);
%WISHART The Wishart matrix.
%     WISHART(N,M,ISREAL) generates an N x N Hermitian Wishart matrix.
%     If ISREAL = 1 then the elements of W are real.
%     If ISREAL = 0 then the elements of W are complex. 
%    
%     N and M are integers while ISREAL equals either 0 or 1.     
%     W is an N x N matrix.
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
%     $Revision: 1.0 $  $Date: 2004/09/10  23:45:18 $

if(isreal==1)
     g = randn(n,m); 
else
     g = (randn(n,m) + i*randn(n,m))/sqrt(2);
end
   
w = (g * g') / m;

