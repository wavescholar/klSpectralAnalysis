function [blurIm,K] = aDiffuse(im,U,S,V,sqrtD,beta0,sizeIm)
% sqrtD, U, S, V are obtained after using binCut0 on
% the affinity matrix. here is a recap:
%
% Acut = A .* ~binCut0; 
% D = sum(Acut, 1)';              % Normalize column sum to one.
% sqrtD = D .^ 0.5;
% Q = (sqrtD .^ -1) * ones(1, length(D));
% Mcut = Q .* Acut .* Q';         % M = D^-0.5 Markov D^0.5
% clear Q
% [U S V] = svd(Mcut);
% S = diag(S);
%
  

  % the normalized symmetric matrix is = U*S*V'
  % to write it as U*S*U' do the following:
  T = diag(V'*U);
  S = S.*T;
  
  % diffusion kernels K : Markov matrix raised to beta0
  K = U*diag(S.^beta0)*U';
  K = diag(sqrtD)*K*diag(sqrtD.^-1);
  % K is a Markov matrix. columns of K sum to 1.
  % so normalization of blurIm by column sums of K is 
  % not necessary.
  blurIm = [sum(K .* (im(:)*ones(1,size(K,2))))]';
  
  % show the blurred image
  figure(713); showIm(reshape(blurIm,sizeIm),[0 255]);
  set(gca,'fontsize',18);
  title(sprintf('Diffused over %d iterations',beta0));
  
  % expect there to be differences when blurIm is compared to
  % the original image, even when beta0 = 1.
  % because A becomes Acut by removing some links
  % so U, S, V correspond to Acut and not A.
  figure(714); showIm(reshape(blurIm-im(:),sizeIm));
  set(gca,'fontsize',18);
  title(sprintf('difference between blurIm and orig image'));
  %imStats(blurIm(:),im(:))
  figure(713);
  
