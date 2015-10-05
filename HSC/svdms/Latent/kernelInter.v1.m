function pixIds = kernelInter(K1,K2,sizeIm,flag)
  
  if nargin < 4
    flag = 0;
  end
  
  % union
  T = K1(:) + K2(:);
  if (flag)
    figure(212);
    subplot(2,2,1);
    showIm(reshape(T,sizeIm));
    title(sprintf('union %d',sum(full(T(:))>0)));
  end
  
  % intersection
  T = K1(:) .* K2(:);
  if (flag)
    subplot(2,2,2);  
    showIm(reshape(T,sizeIm));
    title('intersection');
  end
  
  T(T < 0.01*max(T)) = 0;
  if (flag)
    subplot(2,2,3);  
    showIm(reshape(T,sizeIm));
    %title(sprintf('%2.3f',max(full(T(:)))));
    title('threshold: 0.01*max(T)');
  end
  
  T(T >0) = 1;
  if (flag)
    subplot(2,2,4);  
    showIm(reshape(T,sizeIm));
    title(sprintf('count: %d',sum(full(T(:)))));
  end
  
  pixIds = find(T > 0);
  pixIds = pixIds(:);
  pixIds = pixIds';
