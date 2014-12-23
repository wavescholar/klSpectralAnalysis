function j = manova(n,m1,m2,isreal);
%MANOVA The MANOVA matrix.
%     MANOVA(N,M1,M2,ISREAL) generates an N x N MANOVA matrix.
%     If ISREAL = 1 then the elements of W are real.
%     If ISREAL = 0 then the elements of W are complex. 
%    
%     N , M1 , and M2 are integers while ISREAL equals either 0 or 1.
%     N < M1.
%     J is an N x N matrix.
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

   
w1 = wishart(n,m1,isreal);        % Generate Wishart Matrices       
w2 = wishart(n,m2,isreal);

j =  (w1 + w2) \ w1;

