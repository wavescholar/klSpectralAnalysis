%function [hier,U,S] = svdms(A, nSvd, sizeIm, smallprob, symLock, sym_min, fullimage, indices, protInfo)
function [U, S,hier] = Svdms(A, nSvd, sizeIm)
%
%TEMP: Assumes form where nonzeroes are on the diag of A...
 
if nargin < 3
  sizeIm = [size(A,1), 1];
end

% scaling the affinity matrix
Dorig = full(sum(A,1))';
ok = 0;
if ok
  scl = median(full(Dorig));
  A = A/scl;
end
D = full(sum(A,1))';

hier{1}.A = A;
hier{1}.sizeIm = sizeIm;
%% PROTEIN SPECIFIC ENHANCEMENT
if exist('protInfo')
  hier{1}.protInfo = protInfo;
else
  hier{1}.protInfo = [];
end
hier{1}.newW = speye(size(A));

ok = 0;
if ok
  figure; 
  Z = A; %dz = sum(Z,2); [dzs,dzi] = sort(dz); spy(Z(dzi,dzi)); 
  spy(Z); set(gca,'fontsize',15); xlabel(' '); title(sprintf('Hier: 1')); axis xy;
end

%profile on -detail 'builtin'
lev = 1;
fprintf(2, ' Latent:  lev %d, size %d\n', lev, prod(hier{lev}.sizeIm)); 
pause(0.001);

% put limit on number of iterations at a certain hierarchy size (added by
% vmb, 10-6-10)
oldsize = prod(hier{lev}.sizeIm);
niter = 0; % count number of iterations for each size
maxniter = 1; % maximum number of iterations for a given level size
fctr = 0.0167; % ... added by vmb, 11-12-10: controls lower limit for kernels to
% select (smallest kernel must have stationary prob >= max(kernel
% stationary prob)/fctr (multiply this once by 1.2 to get .02, the starting
% percent limit)




while ( (prod(hier{lev}.sizeIm) > nSvd) && (niter < maxniter)) % [vmb added niter, 10-6-10]
  if lev ==1
    if exist('protInfo')
      logpow =2;        
    else
      logpow =1;      
    end
  else
    logpow = 2;
  end
  lev = lev+1;
  logpow = 2;
  
  if exist('protInfo')
      [hier{lev}.L hier{lev}.A hier{lev}.K hier{lev}.R hier{lev}.st ...
          hier{lev}.W hier{lev}.rbinNhbr hier{lev}.selectId hier{lev}.M ...
          hier{lev}.protInfo] = ...
          buildLatentProtein(hier{lev-1}.A, logpow, hier{lev-1}.sizeIm,hier{lev-1}.protInfo);
      %buildLatent(hier{lev-1}.A, logpow, hier{lev-1}.sizeIm);
      %buildLatentFK(hier{lev-1}.A, logpow, hier{lev-1}.sizeIm);
  else
     
    [hier{lev}.L hier{lev}.A hier{lev}.K hier{lev}.R hier{lev}.st ...
         hier{lev}.W hier{lev}.rbinNhbr hier{lev}.selectId hier{lev}.M] ...
         = buildLatent_orig(hier{lev-1}.A, logpow, hier{lev-1}.sizeIm);
    
  end
  
  
  
  hier{lev}.sizeIm = [size(hier{lev}.A, 1) 1];
  if prod(hier{lev}.sizeIm) == oldsize
      niter = niter+ 1;
      lev = lev-1;
      hier = hier(1:lev);
%       if niter < maxniter/2
%         fctr = fctr*.9;
%       else % aggressively decrease number of kernels allowed!
%           fctr = fctr*.5; 
%       end
    %  disp(['niter for this size: ',num2str(niter),' : ',num2str(maxniter)]);
  else % new size, so start counting again
      niter = 0;
      oldsize = prod(hier{lev}.sizeIm);
  end
  
  ok = 0;
  if ok 
    figure; 
    Z = hier{lev}.A; 
    %dz = sum(Z,2); [dzs,dzi] = sort(dz); spy(Z(dzi,dzi)); 
    spy(Z);
    set(gca,'fontsize',15); xlabel(' '); title(sprintf('Hier: %d',lev)); axis xy;
  end
  
  
  if max(abs(hier{lev}.L(:))) == 0
      disp(['end hierarchy at level ',num2str(lev-1),' bc L = 0s for level ',num2str(lev),', which had ',num2str(prod(hier{lev}.sizeIm)),' kernels.']);
      lev = lev - 1;
      hier = hier(1:lev);
      break;
  end
  
  fprintf(2, ' Latent:  lev %d, size %d\n', lev, prod(hier{lev}.sizeIm));
  pause(0.001);
  
  if hier{lev}.sizeIm(1) < nSvd
    clear hier{lev};
    lev = lev-1;
    %lev
    break;
  end
  
  
end 
nLev = lev;


% % % don't need to compute EVs
 U = []; 
 S = [];
  return; 

%keyboard;
if nLev > 1
  [hier{nLev-1}.U, hier{nLev-1}.S,hier{nLev}.U,hier{nLev}.S,...
   hier{nLev-1}.itr] = coarseFineTPI(hier{nLev}.L,hier{nLev}.A,hier{nLev}.K, ...
                                  nSvd);
  fprintf(2, ' Prolong: lev %d, size %d, nU %d\n', [nLev, size(hier{nLev}.U)]);
  fprintf(2, ' Prolong: lev %d, size %d, nU %d\n', ...
          [nLev-1, size(hier{nLev-1}.U)]);
  pause(0.001);
  for lev = nLev-1:-1:2  
    [hier{lev-1}.U, hier{lev-1}.S,Pp,Qq,hier{lev-1}.itr] ...
        = coarseFineTPI(hier{lev}.L,hier{lev}.A,hier{lev}.K,nSvd, ...
                     hier{lev}.U,hier{lev}.S);
    fprintf(2, ' Prolong: lev %d, size %d, nU %d\n', ...
            [lev-1, size(hier{lev-1}.U)]);
    pause(0.001);
  end
  U = hier{1}.U;
  S = hier{1}.S;
else
   
    D = sum(A,2);
    
    sqrtD    = spdiags(D.^0.5, 0, length(D), length(D));
    sqrtDinv = spdiags(D.^ -0.5, 0, length(D), length(D));
    
    if issparse(A)
      [U,S,V] = svds(sqrtDinv*A*sqrtDinv,...
                     min(nSvd,size(A,2)));
      S = diag(S);
    else
      [U,S,V] = svd(full(sqrtDinv * A * sqrtDinv));
      S = diag(S);
      if nSvd < size(U,2)
        U = U(:, 1:nSvd);
        S = S(1:nSvd);
      end
    end
end
%profile report /var/tmp/summary

return;
