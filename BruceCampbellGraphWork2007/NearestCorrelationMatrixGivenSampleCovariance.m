%%%%%%%%% This code is designed to solve %%%%%%%%%%%%%
%%%%%%%%  min 0.5*<X-G, X-G>
%%%%%%%   s.t. X_ii = b_i>0, i=1,2,...,n
%%%%%%%        X>=tau*I (symmetric and positive semi-definite) %%%%%%%%%%%%%%%
%%%%%%%%
%%%%%%  based on the algorithm  in %%%%%
%%%%%%  ``A Quadratically Convergent Newton Method for %%%
%%%%%%    Computing the Nearest Correlation Matrix %%%%%
%%%%%%%   By Houduo Qi and Defeng Sun  %%%%%%%%%%%%
%%%%%%%   SIAM J. Matrix Anal. Appl. 28 (2006) 360--385.
%%%%%%%  
%%%%%% Last modified date:  September 10, 2009  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% The  input argument is the given symmetric G   %%%%%%%%
%%%%%% The outputs are the optimal primal and dual solutions %%%%%%%%
%%%%%%% Diagonal Preconditioner is added         %%%%%%
%%%%%%% Send your comments and suggestions to    %%%%%%
%%%%%%% hdqi@soton.ac.uk  or matsundf@nus.edu.sg %%%%%%
%%%%%                          %%%%%%%%%%%%%%%
%%%%% Warning: Accuracy may not be guaranteed!!!!! %%%%%%%%

function [X,y] = Correlation_Newton1(G,b,tau)
disp(' ---Newton method starts--- ')
t0 = clock;
[n,m] =size(G);


global  b0

G =(G+G')/2; % make G symmetric
% b0 =ones(n,1);
b0=b;

if nargin==3
   G  =  G - tau*eye(n); % reset G
   b0 =  b0 - tau*ones(n,1); % reset b0   
end

Res_b = zeros(300,1);



y = zeros(n,1);       %Initial point
%y=b0-diag(G);              

Fy = zeros(n,1);

k=0;
f_eval = 0;

Iter_Whole = 200;
Iter_inner = 20; % Maximum number of Line Search in Newton method
maxit = 200; %Maximum number of iterations in PCG
iterk = 0;
Inner = 0;
tol = 1.0e-2; %relative accuracy for CGs

error_tol=1.0e-6; % termination tolerance
sigma_1=1.0e-4; %tolerance in the line search of the Newton method

x0=y;

prec_time = 0;
pcg_time = 0;
eig_time =0;

c = ones(n,1);
%M = diag(c); % Preconditioner to be updated

d =zeros(n,1);

  val_G = sum(sum(G.*G))/2;

 X = G + diag(y);
 X = (X +X')/2;
 
 eig_time0 = clock;
 [P,D] =  eig(X);   %% X= P*diag(D)*P'
 eig_time = eig_time + etime(clock,eig_time0);

% save Cmat.mat C
 %lambda=diag(D);
  P = real(P);
 lambda = real(diag(D));
 if lambda(1) < lambda(n) %the eigenvalues are arranged in the decreasing order
     lambda = lambda(n:-1:1);
    P = P(:,n:-1:1);
 end 
 
 [f0,Fy] = gradient(y,lambda,P,b0,n);
 f = f0;
 f_eval = f_eval + 1;
 b = b0 - Fy;
 norm_b = norm(b);

 Initial_f = val_G - f0;
 fprintf('Newton: Initial Dual Objective Function value==== %d \n', Initial_f)
 fprintf('Newton: Norm of Gradient %d \n',norm_b)
 
 Omega12 = omega_mat(P,lambda,n);
 x0 = y;

 while (norm_b>error_tol & k< Iter_Whole)

  prec_time0 = clock;
   c = precond_matrix(Omega12,P,n); % comment this line for  no preconditioning
  prec_time = prec_time + etime(clock, prec_time0);
  
 pcg_time0 = clock;
 [d,flag,relres,iterk]  =pre_cg(b,tol,maxit,c,Omega12,P,n);
 pcg_time = pcg_time + etime(clock,pcg_time0);
 %d =b0-Fy; gradient direction
 fprintf('Newton: Number of CG Iterations %d \n', iterk)
  
  if (flag~=0); % if CG is unsuccessful
     % d =b0-Fy;
     disp('..... Not a full Newton step......')
  end
 slope = (Fy-b0)'*d; %%% nabla f d
 

    y = x0 + d; %temporary x0+d 
    
      X = G + diag(y);
      X = (X + X')/2; 
     eig_time0 = clock;
     [P,D] = eig(X); % Eig-decomposition: X =P*diag(D)*P^T
     eig_time = eig_time + etime(clock,eig_time0); 
      
    % lambda=diag(D);
      P = real(P);
     lambda = real(diag(D));
     
     if lambda(1) < lambda(n) %the eigenvalues are arranged in the decreasing order
        lambda = lambda(n:-1:1);
        P = P(:,n:-1:1);
     end 
     
     [f,Fy] = gradient(y,lambda,P,b0,n);

     k_inner = 0;
     while(k_inner <=Iter_inner & f> f0 + sigma_1*0.5^k_inner*slope + 1.0e-6)
         k_inner = k_inner+1;
         y = x0 + 0.5^k_inner*d; % backtracking   
         
         X = G + diag(y);
         X = (X + X')/2;
         
         eig_time0 = clock;
         [P,D] = eig(X); % Eig-decomposition: X =P*diag(D)*P^T
          eig_time = eig_time + etime(clock,eig_time0); 
         %lambda = diag(D);
         P = real(P);
         lambda = real(diag(D));
          if lambda(1) < lambda(n) %the eigenvalues are arranged in the decreasing order
           lambda = lambda(n:-1:1);
            P = P(:,n:-1:1);
           end 
         
         
         [f,Fy] = gradient(y,lambda,P,b0,n);
      end % loop for while
      f_eval = f_eval + k_inner+1;
      x0 = y;
      f0 = f;
      
     k=k+1;
     b = b0-Fy;
     norm_b = norm(b);
     fprintf('Newton: Norm of Gradient %d \n',norm_b)

     Res_b(k) = norm_b;
    
     Omega12 = omega_mat(P,lambda,n);

 end %end loop for while i=1;
 

Ip = find(lambda>0);
r = length(Ip);
 
if (r==0)
    X =zeros(n,n);
elseif (r==n)
    X = X;
elseif (r<=n/2)
    lambda1 = lambda(Ip);
    lambda1 = lambda1.^0.5;
    P1 = P(:, 1:r);
    if r >1
        P1 = P1*sparse(diag(lambda1));
        X = P1*P1'; % Optimal solution X* 
    else
        X = lambda1^2*P1*P1';
    end
else 
    
    lambda2 = -lambda(r+1:n);
    lambda2 = lambda2.^0.5;
    P2 = P(:, r+1:n);
    P2 = P2*sparse(diag(lambda2));
    X = X + P2*P2'; % Optimal solution X* 
end
    
 Final_f = val_G-f;
 val_obj = sum(sum((X-G).*(X-G)))/2;
 X = X + tau*eye(n); 
 time_used = etime(clock,t0);
 fprintf('\n')

%fprintf('Newton: Norm of Gradient %d \n',norm_b)
fprintf('Newton: Number of Iterations == %d \n', k)
fprintf('Newton: Number of Function Evaluations == %d \n', f_eval)
fprintf('Newton: Final Dual Objective Function value ========== %d \n',Final_f)
fprintf('Newton: Final Original Objective Function value ====== %d \n', val_obj)
fprintf('Newton: The rank of the Optimal Solution ================= %d \n',r)

fprintf('Newton: computing time for computing preconditioners == %d \n', prec_time)
fprintf('Newton: computing time for linear systems solving (cgs time) ====%d \n', pcg_time)
fprintf('Newton: computing time for  eigenvalue decompostions (calling eig time)==%d \n', eig_time)
fprintf('Newton: computing time used for equal weight calibration ==== =====================%d \n',time_used)



%%% end of the main program



%%%%%%
%%%%%% To generate F(y) %%%%%%%
%%%%%%%

function [f,Fy]= gradient(y,lambda,P,b0,n)
 
%[n,n]=size(P);
f=0.0;
Fy =zeros(n,1);
 
%lambdap=max(0,lambda);
%H =diag(lambdap); %% H =P^T* H^0.5*H^0.5 *P

 P=P';
 i=1;
 while (i<=n)
     P(i,:)=max(lambda(i),0)^0.5*P(i,:);
     i=i+1;
 end
 i=1;
 while (i<=n)
       Fy(i)=P(:,i)'*P(:,i);
 i=i+1;     
 end
 i=1;
 while (i<=n)
     f =f+(max(lambda(i),0))^2;
     i=i+1;
 end
 
f =0.5*f -b0'*y;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% end of gradient.m %%%%%%

%%%%%%%%%%%%%%        %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% To generate the first -order difference of lambda
%%%%%%%

%%%%%%%%%%%%%%        %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% To generate the essential part of the first -order difference of d
%%%%%%%
function Omega12 = omega_mat(P,lambda,n)
%We compute omega only for 1<=|idx|<=n-1
idx.idp = find(lambda>0);
idx.idm = setdiff([1:n],idx.idp);
n =length(lambda);
r = length(idx.idp);
 
if ~isempty(idx.idp)
    if (r == n)
        Omega12 = ones(n,n);
    else
        s = n-r;
        dp = lambda(1:r);
        dn = lambda(r+1:n);
        Omega12 = (dp*ones(1,s))./(abs(dp)*ones(1,s) + ones(r,1)*abs(dn'));
        %  Omega12 = max(1e-15,Omega12);

    end
else
    Omega12 =[];
end

    %%***** perturbation *****
    return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end of omega_mat.m %%%%%%%%%%

%%%%%% PCG method %%%%%%%
%%%%%%% This is exactly the algorithm by  Hestenes and Stiefel (1952)
%%%%%An iterative method to solve A(x) =b  
%%%%%The symmetric positive definite matrix M is a
%%%%%%%%% preconditioner for A. 
%%%%%%  See Pages 527 and 534 of Golub and va Loan (1996)

function [p,flag,relres,iterk] = pre_cg(b,tol,maxit,c,Omega12,P,n);
% Initializations
r = b;  %We take the initial guess x0=0 to save time in calculating A(x0) 
n2b =norm(b);    % norm of b
tolb = tol * n2b;  % relative tolerance 
p = zeros(n,1);
flag=1;
iterk =0;
relres=1000; %%% To give a big value on relres
% Precondition 
z =r./c;  %%%%% z = M\r; here M =diag(c); if M is not the identity matrix 
rz1 = r'*z; 
rz2 = 1; 
d = z;
% CG iteration
for k = 1:maxit
   if k > 1
       beta = rz1/rz2;
       d = z + beta*d;
   end
   %w= Jacobian_matrix(d,Omega,P,n); %w = A(d); 
   w = Jacobian_matrix(d,Omega12,P,n); % W =A(d)
   denom = d'*w;
   iterk =k;
   relres = norm(r)/n2b;              %relative residue = norm(r) / norm(b)
   if denom <= 0 
       sssss=0
       p = d/norm(d); % d is not a descent direction
       break % exit
   else
       alpha = rz1/denom;
       p = p + alpha*d;
       r = r - alpha*w;
   end
   z = r./c; %  z = M\r; here M =diag(c); if M is not the identity matrix ;
   if norm(r) <= tolb % Exit if Hp=b solved within the relative tolerance
       iterk =k;
       relres = norm(r)/n2b;          %relative residue =norm(r) / norm(b)
       flag =0;
       break
   end
   rz2 = rz1;
   rz1 = r'*z;
end

return

%%%%%%%% %%%%%%%%%%%%%%%
%%% end of pre_cg.m%%%%%%%%%%%


%%%%%% To generate the Jacobain product with x: F'(y)(x) %%%%%%%
%%%%%%%

function Ax = Jacobian_matrix(x,Omega12,P,n)

Ax =zeros(n,1);

[r,s] = size(Omega12);
if (r>0)
   % P1 = P(:,1:r);
     H1 = P(:,1:r); 
   if (r< n/2)
         %P2 = P(:,r+1:n);
         %H1=P1;
           i=1;
           while (i<=n)
                H1(i,:) = x(i)*H1(i,:); % H1=diag(x)*P1
                i=i+1;
            end
        
       
        
        Omega12 = Omega12.*(H1'*P(:,r+1:n));


        %H =[(H1'*P1)*P1'+ Omega12*P2';Omega12'*P1']; %%%  H= [Omega o (P^T*diag(x)*P)]*P^T
        H =[(H1'*P(:,1:r))*(P(:,1:r))'+ Omega12*(P(:,r+1:n))';Omega12'*(P(:,1:r))']; %%%  H= [Omega o (P^T*diag(x)*P)]*P^T
        i=1;
        while (i<=n)
            Ax(i)=P(i,:)*H(:,i);
            Ax(i) = Ax(i) + 1.0e-10*x(i); % add a small perturbation
            i=i+1;

        end
    else % r >n/2
        if (r==n)
            Ax =(1+1.0e-10)*x;
        else
            %P2 = P(:,r+1:n);
            H2 =  P(:,r+1:n);
            
             i=1;
            while (i<=n)
                H2(i,:) = x(i)*H2(i,:); % H2=diag(x)*P2
                i=i+1;
            end
            
    
         
            Omega12 = ones(r,s)-Omega12;
            Omega12 = Omega12.*((P(:,1:r))'*H2);

            H =[Omega12* (P(:,r+1:n))';Omega12'*(P(:,1:r))'+ ( (P(:,r+1:n))'*H2)* (P(:,r+1:n))']; %%% Assign H*P' to H= [(ones(n,n)-Omega) o (P^T*diag(x)*P)]*P^T

            i=1;
            while (i<=n)
                Ax(i)= -P(i,:)*H(:,i);
                Ax(i) = x(i) + Ax(i) + 1.0e-10*x(i); % add a small perturbation
                i=i+1;

            end

        end
    end
end
    return
 
%%%%%%%%%%%%%%%
%end of Jacobian_matrix.m%%%

%%%%%% To generate the diagonal preconditioner%%%%%%%
%%%%%%%

function c = precond_matrix(Omega12,P,n)

[r,s] =size(Omega12);
c = ones(n,1);

if (r>0)

    if (r< n/2)
        H = P';
        H = H.*H;

        H12 = H(1:r,:)'*Omega12;
        d =ones(r,1);
        for i=1:n
            c(i) = sum(H(1:r,i))*(d'*H(1:r,i));
            c(i) = c(i) +2.0*(H12(i,:)*H(r+1:n,i));
            if c(i) < 1.0e-8
                c(i) =1.0e-8;
            end
        end
    else  % if r>=n/2, use a complementary formula
        H = P';
        H = H.*H;
        Omega12 = ones(r,s)-Omega12;
        H12 = Omega12*H(r+1:n,:);
        d =ones(s,1);
        dd = ones(n,1);
        
        for i=1:n
            c(i) = sum(H(r+1:n,i))*(d'*H(r+1:n,i));
            c(i) = c(i) + 2.0*(H(1:r,i)'*H12(:,i));
            alpha = sum(H(:,i));
            c(i) = alpha*(H(:,i)'*dd)-c(i);
            if c(i) < 1.0e-8
                c(i) =1.0e-8;
            end
        end

    end
end

return

 
%%%%%%%%%%%%%%%
%end of precond_matrix.m%%%

%



