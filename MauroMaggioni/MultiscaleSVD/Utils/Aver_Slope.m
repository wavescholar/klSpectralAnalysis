function AverSlope = Aver_Slope(X, U)
% finding averages of slope
% Centered difference
% input: 
%   U: a row vector or matrix, denoising each row, size M by N
%   X; space discritization, 1 by N vector
% ouput:
%   AverSlope: average of slopes


[M,N]=size(U); % # of columns

deltaX=X(3:N)-X(1:N-2); % 1 by N-2 vector 
DeltaX=repmat(deltaX, M, 1);  % M by N+1 matrix
Slope=(U(:, 3:N)-U(:,1:N-2))./DeltaX;
AverSlope=sum(Slope,2)/(N-2);

return;