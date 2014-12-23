function [binCut0,cutPar,discComp,U,S,V] = itCut3(Pts,A,binCut0,binNhbr,cutPar)
%function [binCut0,cutPar,discComp,U,S,V] = itCut3(Pts,A,binCut0,binNhbr,cutPar)
%  
% calls: kmeans
%
%
  
FALSE = (0 == 1);
TRUE = ~FALSE;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%  Iterative Cutting. %%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%dLogHalfMin = -0.2; 
%dLogHalfMaxCut = -0.02; 
%maxIts = 1000; 
%maxBasis = size(Pts,1); 

sizeIm         = cutPar.sizeIm;
beta0          = cutPar.beta0 ;
half0          = cutPar.half0 ;
dLogHalfMin    = cutPar.dLogHalfMin ;
dLogHalfMaxCut = cutPar.dLogHalfMaxCut ; 
maxIts         = cutPar.maxIts ; 
maxBasis       = cutPar.maxBasis;
displayON      = cutPar.displayON;
displayStep    = cutPar.displayStep;
useImage       = cutPar.useImage;
svdDone        = cutPar.svdDone;

id0 = 2; 

binCut = zeros(size(binCut0));
%svdDone = FALSE;
%binCut0 = binCut;
voteCut = binCut;
binSupp = binCut;

Acut = A .* ~binCut0; 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%% Cut here for restart after non-convergence. 
if (displayON)
  figure(1); clf; 
end

if ~exist('useImage', 'var') 
  useImage = FALSE; 
end 
 
for its=1:maxIts 
   
  %binCut = zeros(size(binNhbr)); 
  binCut(:) = 0;
  minCut0 = 0; 
  maxCut0 = 0; 
  idCut = 0; 
  change = FALSE; 
  for id = id0:maxBasis 
    %voteCut = zeros(size(binCut)); 
    voteCut(:) = 0;
    
    if ~svdDone 
 
      D = sum(Acut, 1)';              % Normalize column sum to one. 
      sqrtD = D .^ 0.5; 
      Q = (sqrtD .^ -1) * ones(1, length(D)); 
      Mcut = Q .* Acut .* Q';         % M = D^-0.5 Markov D^0.5 
      clear Q 
 
      [U S V] = svd(Mcut); 
      S = diag(S); 
      % S(1:10) 
      S = min(S, 1-10*eps); 
      halfL = -log(2) * log(S.^2).^-1 ; 
      %halfL(2:10) 
      svdDone = TRUE; 
       
    end 
    if (halfL(id) < half0) 
      break; 
    end 
     
    fprintf('id: %d halfL: %f ',id, halfL(id));   
    u = U(:,id); 
    if (displayON & (mod(its,displayStep)==0))    
      figure(1); clf; 
      subplot(2,2,1); 
      if useImage 
	showIm(reshape(u, sizeIm)); 
      else 
	drawEigen(u, Pts); 
      end 
      title(sprintf('U(%d)', id)); 
      subplot(2,2,4); 
      if useImage 
	showIm(reshape(max(binCut0,[],2), sizeIm)); 
      else 
	drawEigen(max(binCut0,[],2), Pts); 
      end 
      title(sprintf('binCut0 (its=%d)', its)); 
      pause(0.0001); 
    end
    
    %tmp = -1/(abs(S(id))*log(abs(S(id)))); 
    tmp = log(2)/(S(id) * log(S(id)) * (beta0 * log(S(id)) - log(2))); 
 
    q = sqrtD.^(-1) .* U(:,id); 
 
    dLogHalf = tmp *(2 * q * q' - ... 
                     S(id)*q.^2 * ones(1, length(q)) - ... 
                     S(id)*(ones(length(q),1) *(q.^2)')); 
 
    % Only consider changing conditional probs for neighbours (but not 
    % auto links, self to self). 
    dLogHalf = dLogHalf .* binNhbr; 
     
    % Non-maximum supression 
    %binSupp = zeros(size(binNhbr)); 
    binSupp(:) = 0;
    
    dLogHalf = reshape(dLogHalf, length(u), length(u)); 
    idx = (dLogHalf(:) < dLogHalfMaxCut) & (~binCut0(:)); 
    if any(idx) 
      tt = find(idx); 
      indIp = ceil(tt/length(u)); 
      indJp = tt - (indIp-1)*length(u); 
      for k=1:length(indIp) 

	suppress = FALSE; 
         
        % Set i0, j0 to be endpoints of edge (indJp(k), indIp(k)) 
        % ordered so that u(i0) > u(j0) 
        if (u(indIp(k)) > u(indJp(k))) 
          i0 = indIp(k); 
          j0 = indJp(k); 
        elseif (u(indJp(k)) > u(indIp(k))) 
          j0 = indIp(k); 
          i0 = indJp(k); 
        else 
          sprintf('Huh? %d %d %f %f\n', indJp(k), indIp(k), u(indIp(k)), dLH) 
          suppress = TRUE; 
        end 
          
        dLH = dLogHalf(j0, i0); 
         
        % Find nhbrs of i0, in current graph with cuts, for which u is 
        % larger than u(i0).  Nonmax suppression occurs if dLogHalf is 
        % still more negative on edges from i0 to these nhbrs. 
        altNhbrs = find(binNhbr(:,i0) & ~binCut0(:, i0) & u(:) > u(i0)); 
        if ~suppress & any(altNhbrs) 
          [mn jmn] = min(dLogHalf(altNhbrs, i0)); jmn = altNhbrs(jmn); 
          if mn < dLH 
            suppress = TRUE; 
          end 
        end 
         
        % Find nhbrs of j0, in current graph with cuts, for which u is 
        % smaller than u(j0).  Nonmax suppression occurs if dLogHalf is 
        % still more negative on edges from j0 to these nhbrs. 
        altNhbrs = find(binNhbr(j0,:) & ~binCut0(j0, :) & (u(j0) > u(:))'); 
        if ~suppress & any(altNhbrs) 
          [mn jmn] = min(dLogHalf(j0, altNhbrs)); jmn = altNhbrs(jmn); 
          if mn < dLH 
            suppress = TRUE; 
          end 
        end 
         
        if suppress 
          binSupp(indJp(k), indIp(k)) = TRUE; 
          binSupp(indIp(k), indJp(k)) = TRUE; 
        end 
      end 
    end  % over idx, all edges with sufficiently small dLogHalf 
     
     
    % Saturate dLogHalf below (for plotting) 
    dLogHalf = max(dLogHalf, dLogHalfMin); 
    dLogHalf = reshape(dLogHalf, length(u), length(u)); 
    [mx ix] = min(dLogHalf,[],2); 
    dQn = diag(dLogHalf(:,ix)); 
     
    if (displayON & (mod(its,displayStep)==0))    
      figure(1);  
      subplot(2,2,2); 
      if useImage 
        showIm(reshape(dQn,sizeIm)); 
      else 
        drawEigen(dQn, Pts); 
      end 
      title(sprintf('dLogHalf(%d)', id)); 
    end
    
    %% Vote for cuts 
 
    % tol = dLogHalfMaxCut; 
    idx = (dLogHalf(:) <= dLogHalfMaxCut) & ~binCut0(:) & ~binSupp(:); 
    if any(idx) 
      tt = find( idx ); 
      indIp = ceil(tt/length(u)); 
      indJp = tt - (indIp-1)*length(u); 
      idx = indIp < indJp; 
      indIp = indIp(idx); 
      indJp = indJp(idx); 
      for k=1:length(indIp) 
	voteCut(indJp(k),indIp(k)) = voteCut(indJp(k),indIp(k)) + 1; 
	voteCut(indIp(k),indJp(k)) = voteCut(indIp(k),indJp(k)) + 1; 
      end 
    end 
    
    if (displayON & (mod(its,displayStep)==0))    
      figure(1); 
      subplot(2,2,3); 
      if useImage 
	showIm(reshape(max(voteCut,[],2), sizeIm)); 
      else 
	drawEigen(max(voteCut,[],2), Pts); 
      end 
      title(sprintf('voteCut(%d)', id)); 
      pause(0.1); 
    end
    
    if any(voteCut(:)>0 & ~binCut0(:)) 
       
      if FALSE 
        minCut = min(dLogHalf(:) .* (voteCut(:)>0) .* ~binCut0(:)) 
        if minCut < minCut0 
          idCut = id; 
          minCut0 = minCut; 
        end 
      else 
        nCut = sum((voteCut(:)>0) .* -dLogHalf(:) .* Acut(:) .* ~binCut0(:));
	fprintf('nCut:  %f ',nCut);   
        if nCut > maxCut0 
          idCut = id; 
          maxCut0 = nCut; 
        end 
      end 
    end 
    fprintf('\n');         
    
  end % over id 
 
  if idCut == 0 
    change = FALSE; 
    %display([change its]); 
    fprintf('\nDONE. its:%d \n',its);
    break; 
  end 
     
  %voteCut = zeros(size(binCut)); 
  voteCut(:) = 0;
  
  id = idCut; 
  fprintf('\nCut selection eigen id: %d halfL %f\n',id, halfL(id));   
  u = U(:,id); 
  
  if (displayON & (mod(its,displayStep)==0))      
    figure(1); 
    subplot(2,2,1); 
    if useImage 
      showIm(reshape(u,sizeIm)); 
    else 
      drawEigen(u, Pts); 
    end 
    title(sprintf('U(%d)', id)); 
    subplot(2,2,4); 
    if useImage 
      showIm(reshape(max(binCut0,[],2), sizeIm)); 
    else 
      drawEigen(max(binCut0,[],2), Pts); 
    end 
    title(sprintf('binCut0 (its=%d)', its)); 
  end
  
  %tmp = -1/(abs(S(id))*log(abs(S(id)))); 
  tmp = log(2)/(S(id) * log(S(id)) * (beta0 * log(S(id)) - log(2))); 
 
  q = sqrtD.^(-1) .* U(:,id); 
 
  dLogHalf = tmp *(2 * q * q' - ... 
                   S(id)*q.^2 * ones(1, length(q)) - ... 
                   S(id)*(ones(length(q),1) *(q.^2)')); 
 
  % Only consider changing conditional probs for neighbours (but not 
  % auto links, self to self). 
  dLogHalf = dLogHalf .* binNhbr; 
   
  % Non-maximum supression 
  %binSupp = zeros(size(binNhbr)); 
  binSupp(:) = 0;
  
  dLogHalf = reshape(dLogHalf, length(u), length(u)); 
  idx = (dLogHalf(:) < dLogHalfMaxCut) & (~binCut0(:)); 
  if any(idx) 
    tt = find(idx); 
    indIp = ceil(tt/length(u)); 
    indJp = tt - (indIp-1)*length(u); 
    for k=1:length(indIp) 
      suppress = FALSE; 
       
      % Set i0, j0 to be endpoints of edge (indJp(k), indIp(k)) 
      % ordered so that u(i0) > u(j0) 
      if (u(indIp(k)) > u(indJp(k))) 
        i0 = indIp(k); 
        j0 = indJp(k); 
      elseif (u(indJp(k)) > u(indIp(k))) 
        j0 = indIp(k); 
        i0 = indJp(k); 
      else 
        sprintf('Huh? %d %d %f %f\n', indJp(k), indIp(k), u(indIp(k)), dLH) 
        suppress = TRUE; 
      end 
       
      dLH = dLogHalf(j0, i0); 
       
      % Find nhbrs of i0, in current graph with cuts, for which u is 
      % larger than u(i0).  Nonmax suppression occurs if dLogHalf is 
      % still more negative on edges from i0 to these nhbrs. 
      altNhbrs = find(binNhbr(:,i0) & ~binCut0(:, i0) & u(:) > u(i0)); 
      if ~suppress & any(altNhbrs) 
        [mn jmn] = min(dLogHalf(altNhbrs, i0)); jmn = altNhbrs(jmn); 
        if mn < dLH 
          suppress = TRUE; 
        end 
      end 
       
      % Find nhbrs of j0, in current graph with cuts, for which u is 
      % smaller than u(j0).  Nonmax suppression occurs if dLogHalf is 
      % still more negative on edges from j0 to these nhbrs. 
      altNhbrs = find(binNhbr(j0,:) & ~binCut0(j0, :) & (u(j0) > u(:))'); 
      if ~suppress & any(altNhbrs) 
        [mn jmn] = min(dLogHalf(j0, altNhbrs)); jmn = altNhbrs(jmn); 
        if mn < dLH 
          suppress = TRUE; 
        end 
      end 
       
      if suppress 
        binSupp(indJp(k), indIp(k)) = TRUE; 
        binSupp(indIp(k), indJp(k)) = TRUE; 
      end 
    end 
  end  % over idx, all edges with sufficiently small dLogHalf 
 
  % Saturate dLogHalf below (for plotting) 
  dLogHalf = max(dLogHalf, dLogHalfMin); 
  dLogHalf = reshape(dLogHalf, length(u), length(u)); 
  [mx ix] = min(dLogHalf,[],2); 
  dQn = diag(dLogHalf(:,ix)); 
  
  if (displayON & (mod(its,displayStep)==0))      
    figure(1); 
    subplot(2,2,2); 
    if useImage 
      showIm(reshape(dQn,sizeIm)); 
    else 
      drawEigen(dQn, Pts); 
    end 
    title(sprintf('dLogHalf(%d)', id)); 
  end
  
  %% Vote for cuts 
 
  % tol = dLogHalfMaxCut; 
  idx = (dLogHalf(:) <= dLogHalfMaxCut) & ~binCut0(:) & ~binSupp(:); 
  if any(idx) 
    tt = find( idx ); 
    indIp = ceil(tt/length(u)); 
    indJp = tt - (indIp-1)*length(u); 
    idx = indIp < indJp; 
    indIp = indIp(idx); 
    indJp = indJp(idx); 
    for k=1:length(indIp) 
      voteCut(indJp(k),indIp(k)) = voteCut(indJp(k),indIp(k)) + 1; 
      voteCut(indIp(k),indJp(k)) = voteCut(indIp(k),indJp(k)) + 1; 
    end 
  end 
  
  if (displayON & (mod(its,displayStep)==0))      
    figure(1); 
    subplot(2,2,3); 
    if useImage 
      showIm(reshape(max(voteCut,[],2), sizeIm)); 
    else 
      drawEigen(max(voteCut,[],2), Pts); 
    end 
    title(sprintf('voteCut')); 
    pause(0.1); 
  end
  
  %% Reduce affinities 
  binCut = (voteCut>0); 
  Acut(binCut) = 0.0; 
   
  change = any(binCut(:) & ~binCut0(:)); 
  fprintf('End of iteration: %d change?: %d\n\n',its,change);
  %display([change its]); 
  if change 
    svdDone = FALSE; 
    % Acut = A .* (~binCut); 
    D = sum(Acut, 1)';        % Normalize column sum to one. 
    sqrtD = D .^ 0.5; 
    Q = (sqrtD .^ -1) * ones(1, length(D)); 
    Mcut = Q .* Acut .* Q';   % M = D^-0.5 Markov D^0.5 
    binCut0 = binCut | binCut0; 
  else 
    break; 
  end 
  
end 

cutPar.its = its;
% Select number of components (see halflife, above) 
dm = find(halfL < 10000); 
%dm = find(halfL < 1000); 
dm = dm(1) - 1 ;
cutPar.dim = dm;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%  END: Iterative Cutting. %%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%  Projection Clustering Using K-means %%%%%%%%%%%%%%%%%%%%% 
%%%%  See Ng, Jordan, and Wiess           %%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%% Params for selecting initial K-means centers 
minRes2 = 0.25;   
MARKOV = FALSE; % MARKOV = TRUE; minRes2 = 0.5; nIts = 64; 
 
 
%%Rerun: rand('state', saveRandState); 
 
n = size(Pts); 
%% Set up data for center selection: 
if MARKOV 
  c = U .* (ones(n(1), 1) * (S .^ nIts)'); 
  E = c; 
  E = ((sum(E.*E, 2).^-0.5) * ones(1, size(E,2))) .* E; 
else 
  E = U(:, 1:dm); 
  E = ((sum(E.*E, 2).^-0.5) * ones(1, dm)) .* E; 
end 
 
%% D0 center selection: 
%% Rerun: rand('state', saveRandState); 
saveRandState = rand('state'); 
 
%% Select initial centers to be nearly orthogonal. 
j = round(0.5 + n(1) * rand(1)); 
j = max(j, 1); j = min(j,n(1)); 
c = E(j,:); 
res = E; 
k = 2; 
while MARKOV | k <= dm 
  res = res - (res * (c(k-1,:)')) * c(k-1,:); 
  nrmRes2 = sum(res .* res, 2); 
  samp = cumsum(nrmRes2 .* (nrmRes2>minRes2)); 
  if samp(n(1)) == 0 | isnan(samp(n(1))) | isinf(samp(n(1))) 
    break; 
  end  
  samp = samp/samp(n(1)); 
  r = rand(1); 
  idx = find(samp>=r); 
  if any(idx) 
    c = [c ; E(idx(1), :)]; 
    k = k+1; 
  else 
    error('Random draw fanned!??!'); 
  end  
end 
k = k-1; 
if k < dm & ~MARKOV 
  fprintf(2,' Got only %d basis elements\n', k); 
  dm = k; 
else 
  fprintf(2,' Got %d basis elements\n', k); 
  dm = k; 
end 
 
%% Call kmeans 
options = foptions; 
[centers options post errlog] = kmeans(c, E, options); 
discComp = post>0; 
 
 
%%%%%%%%%%%%%%%%%%%%%%% SCRATCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
colourSize = 33; 
fcMap = hsv(colourSize);  
leafComp = discComp; 
figure(503); clf; 
if (useImage) 
  qPts = Pts(:, [2,3]);  
  showIm(reshape(Pts(:,1),sizeIm)); 
else 
  qPts = Pts(:, [1,2]); 
  plot(Pts(:,1), Pts(:,2),'.r'); 
end 
hold on; 
for k=1:size(leafComp,2) 
  idx = leafComp(:,k); 
  if any(idx) 
    c = 1+mod(floor(rand*95357),colourSize); 
     
    if FALSE 
    plot(qPts(idx,1), qPts(idx>0,2),'o', 'Color', fcMap(c,:),... 
	 'MarkerFaceColor', fcMap(c,:)); 
    end 
    ht= text(qPts(idx,1), qPts(idx>0,2),sprintf('%d',k), ... 
         'Color', fcMap(c,:), 'FontName', 'Courier', 'FontWeight', 'bold', ... 
             'FontSize', 20); 
	     
 
  else 
    fprintf(2,'empty %d component\n', k); 
  end 
end 
axis equal; axis off; 
hold off; 
 
 
