function j = jsymeig(a)
%JSYMEIG  Jacobian for the symmetric eigenvalue problem
%       J = JSYMEIG(A) returns the Numerical and Theoretical Jacobian for 
%           the symmmetric eigenvalue problem. The numerical answer is 
%           computed by perturbing the symmetric matrix A in the N*(N+1)/2
%           independent directions and computing the change in the
%           eigenvalues and eigenvectos.
%
%     A is a REAL N x N symmetric matrix 
%     Example: Valid Inputs for A are
%              1) G = randn(5); A = (G + G') 
%              2) G = randn(5,10); A = G * G' 
%                                    
%     J is a 1 x 2 row vector
%     J(1) is the theoretical Jacobian computed using the formula in the 
%          third row of the table in Section 12 of Handout # 6
%
%  
%     References:
%     [1] Alan Edelman, Handout 6: Essentials of Finite Random Matrix
%                       Theory, Fall 2004, Course Notes 18.338.
%     [2] Alan Edelman, Random Matrix Eigenvalues.
%
%     Alan Edelman and Raj Rao, Sept. 2004.
%     $Revision: 1.1 $  $Date: 2004/09/28  17:11:18 $

format long
[q,e] = eig(a);  % Compute eigenvalues and eigenvectors
e = diag(e);
epsilon = 1e-7;  % Size of perturbation
n = length(a);   % Size of Matrix
jacmatrix = zeros(n*(n+1)/2);
k = 0; mask = triu( ones(n),1); mask = logical(mask(:));
for i = 1:n,  
    for j = i:n     
        k = k+1;     
        E = zeros(n);     
        E(i,j) = 1; E(j,i) = 1;     
        aa = a + epsilon * E;     
        [qq,ee] = eig(aa);     
        de= (diag(ee)-e)/epsilon;     
        qdq = q'*(qq-q)/epsilon; qdq = qdq(:);     
        jacmatrix(1:n,k) = de;     
        jacmatrix((n+1):end,k) = qdq(mask);  
    end
end

% Numerical answer
j = 1/det(jacmatrix);

% Theoretical answer
[e1,e2] = meshgrid(e,e);
z = abs(e1-e2);
j = abs([j prod( z(mask) ) ]);