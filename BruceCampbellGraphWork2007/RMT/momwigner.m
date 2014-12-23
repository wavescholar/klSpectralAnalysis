function m = momwigner(n,maxmom,isreal);
%MOMWIGNER Moments of the Wigner matrix.
%     M = MOMWIGNER(N,MAXMOM,ISREAL) generates an N x N Wigner matrix.
%                          and computes its first MAXMOM moments.
%     If ISREAL = 1 then the elements of W are real. 
%     If ISREAL = 0 then the elements of W are complex.
%
%     N and MAXMOM are integers while ISREAL equals either 0 or 1.
%     M is a 2 x MAXMOM vector whose first row contains the empirical m
%                moments while its second round contains the empirical
%                moments rounded to the nearest integer.
%        
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
%     Raj Rao, Sept. 2004.
%     $Revision: 1.0 $  $Date: 2004/09/14  09:28:33 $

w = wigner(n,isreal);

e = real(eig(w));

m = [];

for index = 1:maxmom
     m = [m sum((e.^index)/n)];
end

m = [m; round(m)];

  

