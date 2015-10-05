%%=====================================================
%% Run cutDemo, and kill it when it calls itCut
%%
%% random walks across scales
%%  look at: tstLatent.m, mkLatent.m
%%=====================================================
addpath('../MultScale');

SANITY  = FALSE;
svdDone = FALSE;

[sA,D,sL,sM] = normalizeAffty(A);
sqrtD = sparse(diag(D.^0.5));
sqrtDinv = sparse(diag(D .^ -0.5));

[U S V] = svds(sL,50); 
S = diag(S);
S(1:min(length(S),10))

regular = 0;
if ~regular
  %% kernel selection: irregular pyramid
  logpow = 2; % 4;
  sMp = sM;      
  for k = 1:logpow
    sMp = sMp * sMp;
  end
  sLp =  sqrtDinv * sMp * sqrtD ;  % L = D^-0.5 Markov D^0.5
  minRes2 = 0.005;
  %% removing boundary effects while doing kernel selection.
  %% what are the ids along the four walls of the image?
  ids = boundIds(sizeIm,1);
  [selectId,c] = gramSelect(sLp(:,ids),minRes2,1);
  selectId = ids(selectId);
else
  %% kernel selection: regular pyramid
  ok = 1;
  
  if (ok)
    %% sqrt(2) sampling
    selectId = [];
    step = 2;
    for k = 0:step:sizeIm(2)-step+1
      selectId  = [selectId ([1:2:sizeIm(1)] + k*sizeIm(1))];
    end
    for k = 1:step:sizeIm(2)-step+1
      selectId  = [selectId ([2:2:sizeIm(1)] + k*sizeIm(1))];
    end
    selectId = sort(selectId);
  else
    %% (2) sampling
    qq = reshape(1:prod(sizeIm),sizeIm);
    selectId = qq(1:2:size(qq,1),1:2:size(qq,2));
  end
  %pp = zeros(sizeIm);
  %pp(selectId) = 1;
  %figure; showIm(pp)
end

sMd = sM;      
logpow = 2;
%for k = 1:2*logpow-1
for k = 1:logpow-1
  sMd = sMd * sMd;
end
sMp = sMd;

%% Diffusion Kernels
K = sMp(:, selectId(:));

%% mixing coefficient
u0 = D/sum(D);

% W ownership matrix
[st,W] = emKernelFit(u0,K,0); 
% latent space markov transition
R = W'*K; 
figure; plot(sum(R));

% latent space affty
Ar = R * diag(st);
% scaling of latent space affinities.
% not necessary for this expt. however...
Ar = Ar/mean((sum(Ar,1)));

Dr = sum(Ar,2);
[Ur,Sr,Vr] = svd( diag(Dr.^-0.5 ) * Ar * diag(Dr.^-0.5) );
Sr = diag(Sr);
Sr(1:10)

% interpolate, that is extend Ur to fine scale
aU = K*Ur;

% power iteration with gram selection
tt = 20;
%tt = min(50,size(aU,2));
aU = aU(:,1:tt);
for k = 1:2*tt
  aU = sL*aU;
end
[id,aU] = gramFixed(aU,1);
aS = diag(aU'*sL*aU);
[aS,id] = sort(-aS);
aS = -aS;
aU = aU(:,id);

% compare eigenvalues: true vs approximate
figure(113); clf;
plot(aS(1:tt),'xr-'); hold on;
plot(S(1:tt),'x-'); hold on;


plot(Sr(1:size(aU,2)),'xr-'); hold on;
plot(S(1:size(aU,2)).^(2^(logpow-1)),'xb-');
plot(S(1:size(aU,2)).^(2^(logpow)),'xb-');

% and eigenvectors
for k= 1:size(aU,2)
  figure(111); clf;
  subplot(1,2,1);
  showIm(reshape(U(:,k), sizeIm));
  title(sprintf('Orig Basis: %d (%1.4f)',k,S(k)));
  subplot(1,2,2);
  sgn = sign(U(:,k)' * aU(:,k));
  showIm(reshape(sgn*aU(:,k), sizeIm));
  title(sprintf('Approx Cut: %d (%1.4f)',k,aS(k)));  
  pause;
end

%=====================================================
% scale 2

% normalizeAffity(Ar) ;
% sampleKernels;





%=====================================================
%=====================================================
addpath('../MultScale');

SANITY  = FALSE;
svdDone = FALSE;

sA = sparse(A);
D = sum(A, 1)';                 % Normalize column sum to one.
sM = sA * sparse(diag(D.^-1));  % Markov
sqrtD = sparse(diag(D.^0.5));
sqrtDinv = sparse(diag(D .^ -0.5));
sL = sqrtDinv * sA * sqrtDinv;

%D     = sum(A, 1)';              % Normalize column sum to one. 
%sqrtD = D .^ 0.5; 
%Q     = (sqrtD .^ -1) * ones(1, length(D)) ; 
%M     = Q .* A .* Q';         % M = D^-0.5 Markov D^0.5 
%sqrtD    = (diag(D.^0.5));
%sqrtDinv = (diag(D .^ -0.5));
%sM       = M;
%clear Q

% original markov matrix
%[U S V] = svd(M); % M is the laplace matrix not markov
%[U S V] = svd(full(sL)); 
[U S V] = svds(sL); 
S = diag(S);
S(1:min(length(S),10))

%sMp = sM;
%[U,S] = eig(full(sMp));
%S = diag(S);
%[S,idx] = sort(-abs(S));
%S = -S;
%U = U(:,idx);

%=====================================================
%% analysis of the markov matrix
%=====================================================
logpow = 2; % 4;
sMp = sM;      
for k = 1:logpow
  sMp = sMp * sMp;
end
sLp =  sqrtDinv * sMp * sqrtD ;  % L = D^-0.5 Markov D^0.5

% how much diffusion do i see?

%=====================================================
%% kernel selection: irregular pyramid

%% irregular pyramid
%minRes2 = 0.005;
%% removing boundary effects while doing kernel selection.
%% what are the ids along the four walls of the image?
%ids = boundIds(sizeIm,1);
%[selectId,c] = gramSelect(sLp(:,ids),minRes2);
%selectId = ids(selectId);

%=====================================================
%% kernel selection: regular pyramid
pp = zeros(sizeIm);
selectId = [];
step = 2;
for k = 0:step:size(pp,2)-step+1
  selectId  = [selectId ([1:2:size(pp,1)] + k*size(pp,1))];
end
for k = 1:step:size(pp,2)-step+1
  selectId  = [selectId ([2:2:size(pp,1)] + k*size(pp,1))];
end
selectId = sort(selectId);
pp(selectId) = 1;
figure; showIm(pp)

sMd = sM;      
%for k = 1:2*logpow-1
for k = 1:ceil(logpow/2)
  sMd = sMd * sMd;
end
sMp = sMd;

%% Diffusion Kernels
K = sMp(:, selectId);

%% mixing coefficient
u0 = D/sum(D);

[st,W] = emKernelFit(u0,K,1); % W ownership matrix

%% what if i find responsibilites for
%% each individual kernel

%for k = 1:size(sMp,2)
%  fprintf('%d..',k);
%  stK(:,k) = emKernelFit(sMp(:,k),K);
%end

%=====================================================
RR = W'*K;

% Check stationary distribution of R is st.
P = K*st;
Q = K'*diag(P.^-1)*K;
% latent space markov
R = diag(st) * Q;
plot(sum(R))

% latent space affty
Ar = R * diag(st);
%Ar = diag(st) * K' * diag(P.^-1) * K * diag(st) ;

%=====================================================
%% ANALYSIS OF LATENT SPACE AFFINITIES
Dr = sum(Ar,2);
% latent space binNhbr
Ars = Ar;
Tol = 0.01*sum(Ar,2) * ones(1,size(Ar,2));
Ars(Ars < Tol) = 0;% change this M + M'
rbinNhbr = Ars > 0;
rbinNhbr = rbinNhbr - eye(size(rbinNhbr));

% is the affinity matrix diagonal dominant?
% full resolution space
%AA = A;
% latent space
AA = Ar;

aa = diag(AA);
figure; plot(sum(AA > aa*ones(1,size(AA,2))))

AA = diag(st) * K'* diag(P.^-1) * A(:,selectId) ;
MM = AA * diag( sum(AA,2).^-1 ) ;

AA = A(:,selectId)'*A(:,selectId);
MM = AA * diag( sum(AA,2).^-1 ) ;

%=====================================================

%[Ur,Sr] = eig(full(R));
%[Ur,Sr,Vr] = svd(full(R));
[Ur,Sr,Vr] = svd( diag(Dr.^-0.5 ) * Ar * diag(Dr.^-0.5) );

Sr = diag(Sr);
[Sr,idx] = sort(-abs(Sr));
Sr = -Sr;
Ur = Ur(:,idx);

t = length(selectId);
t = min(15,t);
[Sr(1:t) S(1:t).^(2^(logpow-1))]

% look at eigen/rayleigh approximation

% either use the kernels
aU  = K*Ur;

% or the gram-schmidt equivalent
%[selectIdK,cK] = gramSelect(K,minRes2);
%showIm(cK'*cK)
%aU  = cK*Ur;

%for k = 1:length(selectId)
tt = selectId;
for k = 1:length(tt)
  figure(101); clf;
  %subplot(1,2,1);
  %showIm(reshape(cK(:,k), sizeIm));
  %title(sprintf('Gram Basis: %d',k));
  %subplot(1,2,2);
  showIm(reshape(K(:,k), sizeIm));
  title(sprintf('Diffusion Kernel: %d loc:%d',k,tt(k)));  
  pause;
end

% unit vectors
ct  = aU'*aU;
aU  = aU * diag(diag(ct).^-0.5);
% rayleigh on the original markov matrix
aS  = aU'*sM*aU;
daS = diag(aS);
[daS(1:t) S(1:t)]

% approximation after about 10-15 vectors
% does not look so great. this number
% is probably related to the number of kernels used.
for k = 1:length(selectId)
  figure(100); clf;
  subplot(1,2,1);
  showIm(reshape(U(:,k), sizeIm));
  title(sprintf('Orig Basis: %d (%1.4f)',k,S(k)));
  subplot(1,2,2);
  sgn = sign(U(:,k)' * aU(:,k));
  showIm(reshape(sgn*aU(:,k), sizeIm));
  title(sprintf('Approx: %d (%1.4f)',k,daS(k)));  
  pause;
end

% cross talk
pp = U(:,1:length(selectId))'*aU ;
showIm(abs(pp))

%=====================================================
% perturbation analysis, histogram of dLogHalf

dLogHalfMin    = -0.2;            
dLogHalfMaxCut = -0.02;           
beta0 = 40;

id = 2;

% full space vector
dLogHalf = perturbBeta(beta0,S(id),U(:,id),diag(sqrtD)) ;
% fix the diagonal
dLogHalf = dLogHalf .* binNhbr; 
[mx ix]  =  min(dLogHalf,[],2) ; 
dQn      = diag(dLogHalf(:,ix)); 
%showIm(reshape(dQn,sizeIm))
[ht,bn] = histo(dQn);
figure; plot(bn,ht);
figure; plot(diag(A(:,ix)),dQn,'.');
[ht,bn] = histo(A(A>0),64);
figure; plot(bn,ht);

dQn(dQn < dLogHalfMin) = dLogHalfMin;
figure;
showIm(reshape(dQn,sizeIm))

% approximation to full space vector
u = K*Ur(:,id);
u = u/sqrt(u'*u);
s = u'*sM*u;
dLogHalf = perturbBeta(beta0,s,u, diag(sqrtD));

% fix the diagonal
dLogHalf = dLogHalf .* binNhbr; 
[mx ix]  =  min(dLogHalf,[],2) ; 
dQn      = diag(dLogHalf(:,ix)); 
%showIm(reshape(dQn,sizeIm))
[ht,bn] = histo(dQn);
figure; plot(bn,ht);
dQn(dQn < dLogHalfMin) = dLogHalfMin;
figure;
showIm(reshape(dQn,sizeIm))

%=====================================================
% latent space

%dLogHalfr = perturbBeta(beta0,Sr(id),Ur(:,id),sqrt(st)) ;
dLogHalfr = perturbBeta(beta0,Sr(id),Ur(:,id),sqrt(Dr)) ;
dLogHalfr = dLogHalfr .* rbinNhbr; 
[mx ix]   =  min(dLogHalfr,[],2) ; 
dQn       = diag(dLogHalfr(:,ix)); 
dAn       = diag(Ar(:,ix));
[ht,bn] = histo(dQn);
figure; plot(bn,ht);
dQn(dQn < dLogHalfMin) = dLogHalfMin;
figure;
plot(dQn,'o-');
hold on;
x = 1:length(Um(:,id));
y = dLogHalfMaxCut*ones(length(x),1);
plot(x,y,'r-','linewidth',2);


% artificial scaling of latent space affinity matrix
%PP = Ar/max(Ar(:));
PP = Ar/mean(full(sum(Ar,1)));

DD = sum(PP,2);
MM = PP*sparse(diag((DD.^-1)));

figure(112);
subplot(2,2,1); showIm(R);  title('Markov');
subplot(2,2,2); showIm(Ar); title('Affinity');
subplot(2,2,3); showIm(MM); title('Markov approximation');
%subplot(2,2,4); showIm(rbinNhbr); title('Bin Nghbr');
subplot(2,2,4); showIm(MM); title('Markov scaled');

%[Um,Sm] = eig(MM);
[Um,Sm,Vm] = svd(diag(DD.^-0.5 ) * PP * diag(DD.^-0.5));
Sm = diag(Sm);

dLogHalfm = perturbBeta(beta0,Sm(id),Um(:,id),sqrt(DD)) ;
dLogHalfm = dLogHalfm .* rbinNhbr; 
[mx ix]   =  min(dLogHalfm,[],2) ; 
dQn       = diag(dLogHalfm(:,ix)); 

dQn(dQn < dLogHalfMin) = dLogHalfMin;

[ht,bn] = histo(dQn);
figure; plot(bn,ht);
dQn(dQn < dLogHalfMin) = dLogHalfMin;
figure;
plot(dQn,'o-');
hold on;
x = 1:length(Um(:,id));
y = dLogHalfMaxCut*ones(length(x),1);
plot(x,y,'r-','linewidth',2);

%=====================================================
