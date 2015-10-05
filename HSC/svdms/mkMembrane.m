function [im] = mkMembrane(N,fgd)

  if (nargin < 2)
    fgd = 0; % no foreground
  end
  
  seed = round(sum(100*clock));
  rand('seed', seed);
  
  fnameRoot = 'random';
  %sigmar    = 5;
  sigmar    = 3;%1.5;
  
  %%% Build Gaussian filter masks, along with derivatives.
  sigmaSqr  = sigmar*sigmar;
  gFiltSize = 2 * round(3.0 * sigmar) + 1;
  x = [1:gFiltSize] - round((gFiltSize+1)/2);
  gFilt   = exp(- x .* x / (2.0*sigmaSqr));
  gFilt   = gFilt/ sum(gFilt(:));
  
  sz   = 2*N;
  %imo  = rand(sz,sz);
  sigman = 2.5;
  %sigman = 1.5;
  
  imo  = 10 + sigman*randn(sz,sz);
  imor = rconv2sep(imo, gFilt, gFilt);
  imor = imor(1:2:end, 1:2:end);
  %imor= corrDn(imo, gFilt, 'reflect1');
  %imor= corrDn(imor, gFilt', 'reflect1',[2 2]);

  im = imor;

  if (fgd)
    %imi  = rand(sz,sz);
    imi  = 9.5 + sigman*randn(sz,sz);
    imir = rconv2sep(imi, gFilt, gFilt);
    imir = imir(1:2:end,1:2:end);
    %imir= corrDn(imi, gFilt, 'reflect1');
    %imir= corrDn(imi, gFilt', 'reflect1',[2 2]);
    %figure; subplot(1,2,1); showIm(imor); subplot(1,2,2); showIm(imir);
    %figure; showIm(imor+ i*imir); 
    
    %im = imor;
    %%% Build the composite
    %lo = 5;
    %hi = 24;
    %im = imor(lo:hi,lo:hi);
    
    %lo = 3;
    %hi = 14;
    % uncomment if i choose to put a square
    % in the middle
    %lo = 7;
    %hi = 11;
    
    %lo = N/2 -1;
    %hi = lo + N/4;
    
    lo = N/2 - N/4;
    hi = N/2 + N/4;
    
    im(lo:hi,lo:hi) = imir(lo:hi,lo:hi);
  end
  
  %showIm(im)
  
  %pgmWrite(im,'data/random.tmp.pgm');
  %im = pgmRead('data/random.tmp.pgm');
  %imName = [fnameRoot '.' sprintf('%d',seed) '.pgm'];
  
  %figure(219); showIm(im);
