function [U, S] = svdmsProf(A, nSvd, sizeIm)
%function [U, S] = svdmsProf(A, nSvd, sizeIm)
%
%TEMP: Assumes form where nonzeroes are on the diag of A...
 
if nargin < 3
  sizeIm = [size(A,1), 1];
end

% scaling the affinity matrix
Dorig = full(sum(A,1))';
scl = median(full(Dorig));
A = A/scl;
D = full(sum(A,1))';

hier{1}.A = A;
hier{1}.sizeIm = sizeIm;

%profile on -detail 'builtin'
lev = 1;
fprintf(2, ' Latent:  lev %d, size %d\n', lev, prod(hier{lev}.sizeIm)); 
pause(0.001);
while prod(hier{lev}.sizeIm) > 180
  if lev ==1
    logpow =1;
  else
    logpow = 2;
  end
  lev = lev+1;
  [hier{lev}.L hier{lev}.A hier{lev}.K hier{lev}.R hier{lev}.st ...
   hier{lev}.W hier{lev}.rbinNhbr sId1 sMp1] = ...
      buildLatent(hier{lev-1}.A, logpow, hier{lev-1}.sizeIm);
      %buildLatentFK(hier{lev-1}.A, logpow, hier{lev-1}.sizeIm);
  hier{lev}.sizeIm = [size(hier{lev}.A, 1) 1];
  fprintf(2, ' Latent:  lev %d, size %d\n', lev, prod(hier{lev}.sizeIm));
  pause(0.001);
  if hier{lev}.sizeIm(1) < nSvd
    clear hier{lev};
    lev = lev-1;
    break;
  end
end 
nLev = lev;

if nLev > 1
  [hier{nLev-1}.U, hier{nLev-1}.S,hier{nLev}.U,hier{nLev}.S,...
   hier{nLev-1}.itr] = coarseFine(hier{nLev}.L,hier{nLev}.A,hier{nLev}.K, ...
                                  nSvd);
  fprintf(2, ' Prolong: lev %d, size %d, nU %d\n', [nLev, size(hier{nLev}.U)]);
  fprintf(2, ' Prolong: lev %d, size %d, nU %d\n', ...
          [nLev-1, size(hier{nLev-1}.U)]);
  pause(0.001);
  for lev = nLev-1:-1:2  
    [hier{lev-1}.U, hier{lev-1}.S,Pp,Qq,hier{lev-1}.itr] ...
        = coarseFine(hier{lev}.L,hier{lev}.A,hier{lev}.K,nSvd, ...
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

fName = sprintf('/var/tmp/summary%d',sizeIm(1));

%profile report /var/tmp/summary
%profile report 'fName'

return;
