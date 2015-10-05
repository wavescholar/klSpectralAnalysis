addpath('Latent');
addpath('newdata');

clear all
close all

FALSE = (0 == 1);
TRUE = ~FALSE;

MAKE_AFF = FALSE;

Nid = [
    16 %1
    24 %2 
    32 %3
    48 %4
    63 %5
    64 %6
    65 %7
    80 %8
    100 %9
    120 %10
    124 %11
    127 %12
    128 %13
    129 %14
    130 %15
    160 %16
    255 %17
    256 %18
    257 %19
    300 %20
    362 %21
    511 %22
    512 %23
    513 %24
];

%Nid(1)  = 16;
%Nid(2)  = 24;
%Nid(3)  = 32;
%Nid(4)  = 48;
%Nid(5)  = 64;
%Nid(6)  = 80;
%Nid(7)  = 100; 
%Nid(8)  = 120;
%Nid(9)  = 124;
%Nid(10) = 128;
%Nid(11) = 130;
%Nid(12) = 160;
%Nid(13) = 256;
%Nid(14) = 300;  % roughly 360*240 pixels
%Nid(15) = 362;  % roughly (256^2)*2 pixels (Laptop OUT OF MEMORY)
%Nid(16) = 512;  

Tres = zeros(length(Nid),11);
%Tres(:,1) = Nid(:);

for cid= 16
  % cid=15;
  %cid=4;

  N = Nid(cid);
  fprintf('\nN: %d \n',N);
  Tres(cid,1) = N;
  
  %dLogHalfAll = [];
  
  %=== affty matrix parameters
  afftyPar.sizeIm  = [N N];
  afftyPar.dsThres = 1.1;
  afftyPar.dsSupp  = 3.1; 
  afftyPar.rho     = 1.5; 
  beta0 = 40;
  half0 = beta0/4;  
  id0 = 2;

  for kk = 1:10,

    fprintf('%d_%d \n',N,kk);
    %=== create membrane. load from cache
    if TRUE
      %load(sprintf('im_%d',N));
      load(sprintf('im_%d_%d',N,kk));      
    else
      % or create a new one
      im = mkMembrane(N,0); % no foreground % save(sprintf('im_%d',N),'im');
      MAKE_AFF = TRUE;
      %im = mkMembrane(N,1); % with foreground
      
      % show the membrance
      %figure(2); clf; showIm(im); pause(.1);
      %save(sprintf('data/im_%d',N),'im');  
      save(sprintf('im_%d_%d',N,kk),'im');        
    end
    
    sizeIm = [N N];%size(im);
    
    %Pts = ones(prod(sizeIm),2);
    %Pts(:,1) = im(:);
    %figure(1);clf;showIm(im);
    
    % affinity matrix. 
    % load from the cache
    if ~MAKE_AFF
      %load(sprintf('data/affty_%d',N));
      load(sprintf('affty_%d_%d',N,kk));      
    else
      % or create a new affinity matrix
      %profile on -detail 'builtin'
      %tic
      A = shiftAffty(im,afftyPar.rho); 
      %toc
      %save(sprintf('data/affty_%d',N),'A');
      save(sprintf('affty_%d_%d',N,kk),'A');      
    end
    
    nSVD = 51;
    tic; 
    [U,S] = svdms(A,nSVD,sizeIm); 
    Tres(cid,kk+1) = toc; 
    fprintf('=> sec: %2.2f \n',Tres(cid,kk+1));   
  end % kk

  % time stats
  minT = min(Tres(cid,2:end));
  maxT = max(Tres(cid,2:end));
  medT = median(Tres(cid,2:end));
  avgT = mean(Tres(cid,2:end));
  stdT = std(Tres(cid,2:end));
  fprintf('N: %d min: %2.2f max: %2.2f median: %2.2f avg: %2.2f std: %2.2f \n',...
	  N,minT,maxT,medT,avgT,stdT);
  
  %for kk = 1:10
  %  fprintf('%2.2f ',Tres(cid,kk+1));
  %end
  
  fprintf('\n');
end

ok = 0;
if ok
figure(101); clf; 
%if ~cid
plot(Tres(1:cid,1).^2, Tres(1:cid,2), '-*b','linewidth',2);
grid on;
set(gca,'fontsize',15);
xlabel('Nodes');
ylabel('Seconds');
title('Time for the leading 51 eigenpairs');
axis tight;
end
