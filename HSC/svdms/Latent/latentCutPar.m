function rcutPar = latentCutPar(maxBasis,sizeIm)   


  FALSE = (0 == 1);
  TRUE = ~FALSE;
  
  
  rcutPar.beta0 = 30;%30;%30.0;  % 10 worked for others
  
  rcutPar.half0          = rcutPar.beta0/10; %4; % /10  or half0 = 10;
  rcutPar.dLogHalfMin    = -0.2;             % dont change this
  rcutPar.dLogHalfMaxCut = -0.1;%-0.01;%-0.07;     % dont change this
  
  %cutPar.dLogHalfMaxCut = -0.02*2^(2*logpow-1) % dont change this
  
  rcutPar.maxIts   = 1;        % dont change this
  rcutPar.maxBasis = maxBasis; %size(Pts,1); % fixed
  rcutPar.its      = 0; % set by itCut3
  rcutPar.dim      = 0; % set by itCut3
  
  rcutPar.sizeIm      = sizeIm;%[size(Ar,1) 1];
  rcutPar.svdDone     = FALSE;
  rcutPar.nSVD        = 10;    % Number of svd vectors to compute first
  rcutPar.useImage    = TRUE;
  rcutPar.displayON   = 0; % set this to 0 for no display
  rcutPar.displayStep = 1; % which iterations to show
  rcutPar.optmethod = 'svd';
  rcutPar.idCut     = [];
  
  rcutPar.adapt     = FALSE;

  
  % works well for
  % logpow = 2
  % beta0 = 30
  % half0 = beta0/10
  % dLogHalfMaxCut = -0.1
