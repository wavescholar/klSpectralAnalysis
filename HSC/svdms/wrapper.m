addpath('Latent');
addpath('data');

clear all
close all

FALSE = (0 == 1);
TRUE = ~FALSE;

MAKE_AFF = FALSE;

Nid(1)  = 16;
Nid(2)  = 24;
Nid(3)  = 32;
Nid(4)  = 48;
Nid(5)  = 64;
Nid(6)  = 80;
Nid(7)  = 100; 
Nid(8)  = 120;
Nid(9)  = 124;
Nid(10) = 128;
Nid(11) = 130;
Nid(12) = 160;
Nid(13) = 256;
Nid(14) = 300;  % roughly 360*240 pixels
Nid(15) = 362;  % roughly (256^2)*2 pixels (Laptop OUT OF MEMORY)
Nid(16) = 512;  

if TRUE %FALSE
  Tres = zeros(length(Nid),3);
  Tres(:,1) = Nid(:);
else
  load 'Tres'
  figure(1); clf; 
  % figure(1); hold on;
  plot(Tres(1:14,1).^2, Tres(1:14,2), '-*b');
end
T = zeros(length(Nid),2);
T(:,1) = Nid(:);


% for cid=1:15
% cid=15;
cid=4;

N = Nid(cid)
Tres(cid,1) = N;

dLogHalfAll = [];

%=== affty matrix parameters
afftyPar.sizeIm  = [N N];
afftyPar.dsThres = 1.1;
afftyPar.dsSupp  = 3.1; 
afftyPar.rho     = 1.5; 
beta0 = 40;
half0 = beta0/4;  
id0 = 2;

%=== create membrane. load from cache
if FALSE
  load(sprintf('im_%d',N));
else
  % or create a new one
  im = mkMembrane(N,0); % no foreground % save(sprintf('im_%d',N),'im');
  MAKE_AFF = TRUE;
  %im = mkMembrane(N,1); % with foreground

  % show the membrance
  %figure(2); clf; showIm(im); pause(.1);
end

sizeIm = [N N];%size(im);

%Pts = ones(prod(sizeIm),2);
%Pts(:,1) = im(:);
%figure(1);clf;showIm(im);

% affinity matrix. 
% load from the cache
if ~MAKE_AFF
  load(sprintf('affty_%d',N));
else
  % or create a new affinity matrix
  %profile on -detail 'builtin'
  %tic
  A = shiftAffty(im,afftyPar.rho); 
  %toc
  save(sprintf('../../svdms/data/affty_%d',N),'A');
end

nSVD = 51;
tic; [U,S] = svdms(A,nSVD,sizeIm); toc
