function C=validcorr(A, niters)
%Based on Nicholas J. Higham, 2002, 
%"Computing the Nearest Correlation Matrix - A Problem from Finance"
%http://eprints.ma.man.ac.uk/232/01/covered/MIMS_ep2006_70.pdf
S=zeros(size(A));
Y=A;
for k=1:length(A)*niters
    R=Y-S;
    [Q D]=eig(R);
    D
    X=Q*max(D,0)*Q';
    S=X-R;
    Y=X-diag(diag(X))+eye(length(X));
end
C=X-diag(diag(X))+eye(length(X));
C