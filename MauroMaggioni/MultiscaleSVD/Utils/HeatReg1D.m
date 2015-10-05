function U=HeatReg1D(X, U, lambda, n)

% Tychonoff regularization using grad^2
% for 1-D with Neumann condition
% Gauss-Jacobi iterative method
% input: 
%   U: a row vector or matrix, denoising each row, size M by N
%   X; space discritization, 1 by N vector
%   n: the number of iteration
%   lambda: weight on the fidelity term
% ouput: 
%    U_new: denoised row vector or matrix, size M by N
% Date: Mar 5, 2008


[M,N]=size(U); 

%Preprocessing: removing Inf and NaN
badPts=find((isinf(U))|(isnan(U)));
U(badPts)=0;
Lambda=lambda*ones(M,N)*max(max(abs(U)));
Lambda(badPts)=0;

U0=U;
Unew=U;
deltaX=X(2:N)-X(1:N-1); % 1 by N-1 vector 
if min(deltaX)<1e-8,
    return;
end;
deltaX=1./[deltaX(1), deltaX, deltaX(N-1)]; % Neumann extension, 1 by N+1 vector
DeltaX=repmat(deltaX, M, 1);                % M by N+1 matrix                                
Divisor = (DeltaX(:, 1:N)+DeltaX(:, 2:N+1)+Lambda);
LambdaU0 = Lambda.*U0./Divisor;
DeltaX1 = DeltaX(:,1:N)./Divisor;
DeltaX2 = DeltaX(:, 2:N+1)./Divisor;
clear DeltaX;

for i=1:n,
    %Uext=[U(:,1), U, U(:, N)]; % Neumann extension, M by N+2 matrix           
    Uext1 = [U(:,1) U(:,1:N-1)];
    Uext2 = [U(:,2:N ), U(:,N)];
    U = Uext1.*DeltaX1+Uext2.*DeltaX2+LambdaU0;
end

return;