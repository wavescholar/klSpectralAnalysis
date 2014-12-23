%jacobian2by2.m
%Code 8.1 of Random Eigenvalues by Alan Edelman

%Experiment:    Compute the Jacobian of a 2x2 matrix function
%Comment:       Symbolic tools are not perfect.  The author
%               exercised care in choosing  the variables.

syms p q r s a b c d t e1 e2
X=[p q ; r s]; A=[a b;c d];

%% Compute Jacobians

Y=X^2;                J=jacobian(Y(:),X(:)), JAC_square  =factor(det(J)) 
Y=X^3;                J=jacobian(Y(:),X(:)), JAC_cube    =factor(det(J))
Y=inv(X);             J=jacobian(Y(:),X(:)), JAC_inv     =factor(det(J)) 
Y=A*X;                J=jacobian(Y(:),X(:)), JAC_linear  =factor(det(J))
Y=[p q;r/p det(X)/p]; J=jacobian(Y(:),X(:)), JAC_lu      =factor(det(J))

x=[p s r];y=[sqrt(p) sqrt(s) r/(sqrt(p)*sqrt(s))];          
                      J=jacobian(y,x),       JAC_DMD     =factor(det(J))
                      
x=[p s]; y=[ sqrt(p^2+s^2) atan(s/p)];                      
                      J=jacobian(y,x),       JAC_notrace =factor(det(J))
                      
Q=[cos(t) -sin(t); sin(t) cos(t)];
D=[e1 0;0 e2];Y=Q*D*Q.';
y=[Y(1,1) Y(2,2) Y(1,2)]; x=[t e1 e2];                      
                      J=jacobian(y,x),       JAC_symeig  =simplify(det(J))
X=[p s;s r]; Y=A.'*X*A; 
y=[Y(1,1) Y(2,2) Y(1,2)]; x=[p r s];
                      J=jacobian(y,x),       JAC_symcong =factor(det(J))