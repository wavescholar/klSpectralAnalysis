function [hier,U,S] = svdms_msm(A, nSvd, nrounds)
%function [U, S,hier] = svdms(A, nSvd, nrounds)
% nSvd: max number of clusters to obtain
% nrounds: max number of hierarchy rounds

%
%TEMP: Assumes form where nonzeroes are on the diag of A...
        %BBCREVISIT This needs explanation
if nargin < 3
  nrounds = 100000000; 
end

sizeIm = [size(A,1), 1];

% scaling the affinity matrix
Dorig = full(sum(A,1))';
ok = 0;
if ok
  scl = median(full(Dorig));                 %BBCREVISIT This needs explanation, and if useful add this to the parameter set

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

hier{1}.selectId = [1:sizeIm]; 
hier{1}.SelectId = [1:sizeIm]; 

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

smallprob = 0;

while ( (prod(hier{lev}.sizeIm) > nSvd) && (niter < maxniter)) % [vmb added niter, 10-6-10]
  
  lev = lev+1;
  if lev == 2
      logpow = 1; % only do sM^2 for largest matrix
  else
      logpow = 2;
  end
  
  if lev > nrounds
      break;
  end
  
  if exist('protInfo')
      [hier{lev}.L hier{lev}.A hier{lev}.K hier{lev}.R hier{lev}.st ...
          hier{lev}.W hier{lev}.rbinNhbr hier{lev}.selectId hier{lev}.M ...
          hier{lev}.protInfo] = ...
          buildLatentProtein(hier{lev-1}.A, logpow, hier{lev-1}.sizeIm,hier{lev-1}.protInfo);

  else
      fctr = min(fctr*1.2,.05)% min(fctr*1.2,.55); %min(fctr + .05, .55); %max(fctr - 10,10);
      kfctr = .3; % remove all kernels who do not have a member that belongs to them by > kfctr percent

    [hier{lev}.L hier{lev}.A hier{lev}.K hier{lev}.R hier{lev}.st ...
         hier{lev}.W hier{lev}.rbinNhbr hier{lev}.selectId hier{lev}.M] ...
         = buildLatent_orig_msm2(hier{lev-1}.A);
     
     if lev > 1
         hier{lev}.newW = hier{lev-1}.newW*hier{lev}.W;
         [tr cl] = max(hier{lev}.newW,[],2);
         hier{lev}.cl = cl;
         
         ncl = max(hier{lev}.cl);
         clsizes = zeros(ncl,1);
         for k = 1:ncl
             clsizes(k) = size(find(hier{lev}.cl == k),1);
         end
         [a,sortclsizes] = sort(clsizes,'descend');
         sortclsizes = sortclsizes(1:size(find(a),1));
         
         hier{lev}.A = hier{lev}.A(sortclsizes,sortclsizes);
         
         hier{lev}.R = hier{lev}.R(sortclsizes,sortclsizes);
         hier{lev}.st = hier{lev}.st(sortclsizes);
         hier{lev}.selectId = hier{lev}.selectId(sortclsizes);
         hier{lev}.M = hier{lev}.M(sortclsizes,sortclsizes);
          hier{lev}.newW = hier{lev}.newW(:,sortclsizes);
         
         ncl = size(sortclsizes,1);
         sortcl = zeros(size(hier{lev}.cl));
         for k = 1:ncl
             thisk = sortclsizes(k);
             u = find(hier{lev}.cl == thisk);
             sortcl(u) = k;
         end
         
         hier{lev}.cl = sortcl;
         hier{lev}.clsizes = a;
         
         hier{lev}.sortindices = sortclsizes;
     end
     
      
     if lev > 2
        hier{lev}.SelectId = hier{lev-1}.SelectId(hier{lev}.selectId);
     else
         hier{lev}.SelectId = hier{lev}.selectId;
     end
     
  end
  
  
  
  hier{lev}.sizeIm = [size(hier{lev}.A, 1) 1];
  if prod(hier{lev}.sizeIm) == oldsize
      niter = niter+ 1;
      lev = lev-1;
      hier = hier(1:lev);
      
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

return;
