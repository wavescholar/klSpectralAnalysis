function y=permanent(A)

%PERMANENT(A) Calculate the permanent of a square matrix A, which is
%defined as the analog of determinant where the signs of each term in 
%summation was removed.
%Example, the matrix 
%A=[1 2
%   3 4]
% is p(A)=1x4+2x3=10.

%written by C.Xu, Nov.,2008, Hangzhou,China. All rights reserved.  


[m,n]=size(A);
%Restrict A to be square
if (m~=n)
    error('A must be square');
end
   if n==1
       y=A;
   else 
       for k=2:n
    P=ones(1,k);
   
    for i=1:k
        SubA=A([1:k-1],[1:i-1 i+1:k]);
        P(i)=permanent(SubA);
    end
    y=A(k,1:k)*P';
       end
   end
   