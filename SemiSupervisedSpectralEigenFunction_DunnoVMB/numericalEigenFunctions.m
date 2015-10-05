function [bins,g,lambdas,pp]=numericalEigenFunctions(DATA,SIGMA,EPSILON)
  %
  %  [bins,g,lambdas,pp]=numericalEigenFunctions(DATA,SIGMA,EPSILON)
  %
  % Function that computes approximate numerical values of eigenfunctions
  % of 1D density at a set of discrete locations. This function solves
  % for g and \sigma in eqn. 2 on page 4 of the NIPS 2009 paper:  
  % "Semi-supervised learning in gigantic image collections" 
  % by R. Fergus, Y. Weiss and A. Torralba. 
  %
  % Inputs:
  %   1. DATA - nPoints x 1 real vector of data from a single dimension.
  %   2. SIGMA - 1 x 1 real -- scaling term in affinity  
  %   3. EPSILON (optional) - 1 x 1 real -- constant added to histogram
  %   to avoid zero density. Default value = 0.1. 
  %
  % Outputs:
  %   1. bins - NUM_BINS x 1 vector -- set of discrete locations at which
  %   the eigenfunctions are computed.
  %   2. g - NUM_BINS x NUM_BINS matrix of (approximate numerical)
  %   eigenfunctions. Each column is a different
  %   eigenfunction. Values of each eigenfunction are computed at the
  %   locations bins (output 1). Equivalent to g in eqn. of paper.
  %   3. lambdas - NUM_BINS x 1 vector of eigenvalues for each
  %   eigenfunction. Equivalent to \sigma in eqn. 2 of paper.
  %   4. pp - NUM_BINS x 1 vector -- histogram of data, with bins
  %   centered at bins (output 1). 
  %
  %    
  % This function should be called from eigenfunctions.m    
  % 
  % Version 1.0. Rob Fergus & Yair Weiss. 11/30/09.  
  
    
  
  % Default constant added to histogram to avoid regions of zero
  % density. See first (non-bulleted) paragraph of p.5 of the paper.
  if nargin==2
    EPSILON=0.1; % do not make too small otherwise eigenfunctions will
                 % become unstable.
  end

  % Maximum number of histogram bins to use. 
  % Solution is not that sensitive to this.
  NUM_BINS = 50;
  
  % measure eigenvalues of eigenfunctions on original density (=1) or the
  % version with EPSILON added (=0). =0 option is less accurate but more stable.
  EIGENVALUES_ORIGINAL_PDF = 0;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Vectorize data (it should be 1D anyway)
  DATA=DATA(:);
  
  % Set #bins: if we have few points then use less than NUM_BINS bins,  
  % otherwise use NUM_BINS.
  nBins=length(DATA)/10;
  nBins=min(nBins,NUM_BINS);
  
  % Now build histogram 
  [pp,bins]=hist(DATA,nBins);
  bins=bins(:);
  pp=pp/sum(pp); % normalize
  pOld=pp; % store for future use (before constant is added)
  pp=pp+EPSILON; % avoid zeros in pdf by adding constant
  pp=pp/sum(pp); % normalize again
    
  % build \tilde{W} in eqn. 2 of paper (page 4). 
  D=dist_mat([bins bins],[bins bins])/sqrt(2);  % just a hack to work on 1D vectors
  W=exp(-0.5*D.^2/SIGMA^2); % W is \tilde{W} in eqn.2.
  
  % P in eqn. 2
  P=diag(pp);
  
  % hatD is \hat{D} in eqn.2 
  hatD=diag(sum(P*W,1));
  
  % PWP is P*\tilde{W}*P in eqn. 2 
  % tildeD is \tilde{D} in eqn. 2 
  PWP=P*W*P;
  tildeD=diag(sum(PWP,1));
  
  % form left-hand-side of eqn. 2
  L=tildeD-PWP;
  
  % we want to minimize x^T L x / x^T (P*hatD) x
  IP=inv(sqrt(P*hatD));
  L2=IP*L*IP;

  %%% Solve for g's in equ. 2
  [uu,ss,vv]=svd(L2);
  g=IP*uu; 
  % g's are numerical approximation of eigenfunctions
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Now measure eigenvalues of the eigenfunctions g
  
  if EIGENVALUES_ORIGINAL_PDF
    P=diag(pOld); % measure with original pdf
    % this is the version before the EPSILON constant was added
    % it should be more accurate using this, but can also be a bit more unstable.

    % recompute terms in eqn.2
    hatD=diag(sum(P*W,1));
    PWP=P*W*P;
    tildeD=diag(sum(PWP,1));
    L=tildeD-PWP;
  end

  % measure eigenvalues
  lambdas=diag(g'*L*g);
  lambdas=lambdas./diag(g'*P*hatD*g);

