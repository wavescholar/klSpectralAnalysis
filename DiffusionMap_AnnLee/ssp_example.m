% Example code for diffusion maps (diffuse.m) 
% and diffusion coarse-graining (diffusion_kmeans.m)
% using real model SSP spectra from Bruzual & Charlot (2003)
%
% Ann Lee, March 2009. Last modified, JWR: 3/23/2009

%--------------------------------------------------------------------
% LOAD DATA
%--------------------------------------------------------------------
close('all');  % close all figures
disp('Loading SSP data')
ssps = load('ssps.txt'); % load SSP spectra
params = load('ssp_params.txt'); % load SSP age/Z
t = params(:,1); % SSP ages
Z = params(:,2); % SSP Zs
lam = load('ssps_wl.txt'); % load wavelength dispersion
[n,p]=size(ssps); % Data(n,p), where n=#observations, p=#variables
lam0 = find(lam==4020); % normalization wavelength
norm = ssps(:,lam0); % keep track of normalization const.
ssps = ssps./repmat(norm,1,p); % normalize all SSPs at lam0

D=squareform(pdist(ssps)); % pairwise distances, n-by-n matrix

%--------------------------------------------------------------------
% SET PARAMETERS IN MODEL
%--------------------------------------------------------------------
eps_val=500/4;
flag_t=0; %flag_t=0 => Default: multi-scale geometry
if flag_t
    t=3;  % fixed time scale  
end

%--------------------------------------------------------------------
% EIGENDECOMPOSITION
%--------------------------------------------------------------------
% Call function "diffuse.m"
disp('Computing diffuson coefficients')
[X, eigenvals, psi, phi] = diffuse(D,eps_val);


%--------------------------------------------------------------------
% DIFFUSION K-MEANS
%--------------------------------------------------------------------
disp('Running diffusion K-means')
k=45;  % number of clusters
Niter=100;
phi0=phi(:,1);
[idx, C, ERR, DX] = diffusion_kmeans(X, k, phi0, Niter);
for jj=1:k, % define prototype parameters
  protospecnorm(jj,:) =  sum(repmat(phi0(idx==jj,1)/sum(phi0(idx==jj,1)),1,p).*ssps(idx==jj,:),1);
  protoages(jj,:) = sum(phi0(idx==jj,1)/sum(phi0(idx==jj,1)).*t(idx==jj),1);
  protometals(jj,:) = sum(phi0(idx==jj,1)/sum(phi0(idx==jj,1)).*Z(idx==jj),1);
  norms(jj,:) = sum(phi0(idx==jj,1)/sum(phi0(idx==jj,1)).*norm(idx==jj),1);
end
protospec= protospecnorm.*repmat(norms,1,p);

%--------------------------------------------------------------------
% PLOT RESULTS
%--------------------------------------------------------------------
figure, % 3-dimensional diffusion map of 1278 SSPs
scatter3(X(:,1),X(:,2),X(:,3),20,idx,'filled'); 
title('Embedding of SSPs with first 3 diffusion coordinates plus K-means centroids');
xlabel('X_1'); ylabel('X_2'); zlabel('X_3');
hold on
scatter3(C(:,1),C(:,2),C(:,3),45,'k','filled'); 


% Cid Fernandes et al.(2005) parameters
tCF = [0.001 0.00316 0.00501 0.01 0.02512 0.04 0.10152 0.28612 0.64054 ...
       0.90479 1.434 2.5 5 11 13]*10^9;
ZCF = [0.0040  0.0200  0.0500];
idxCF=[]; % find SSP spectra used by Cid Fernandes et al.(2005)
for ii=1:15,
  for jj=1:3,
    idxCF = [idxCF find(t==tCF(ii) & Z==ZCF(jj))];
  end
end

% plot normalized prototype spectra, colored by log age
cmap = colormap;
CFcolor = ceil((log10(t(idxCF))-min(log10(protoages)))/(max(log10(protoages))-min(log10(protoages)))*63)+1;
dmapcolor = ceil((log10(protoages)-min(log10(protoages)))/(max(log10(protoages))-min(log10(protoages)))*63)+1;

figure, % prototype spectra, colored by log age
subplot(2,1,1)
for ii=1:45,
  plot(lam,ssps(idxCF(ii),:),'Color',cmap(CFcolor(ii),:),'LineWidth',.01)
  hold on
end
xlabel('lambda'); ylabel('Normalized flux');ylim([0 3]);
title('Normalized Cid Fernandes et al.(2005) prototype spectra, K=45, colored by log age');
subplot(2,1,2)
for ii=1:k,
  plot(lam,protospecnorm(ii,:),'Color',cmap(dmapcolor(ii),:),'LineWidth',.01)
  hold on
end
xlabel('lambda'); ylabel('Normalized flux'); ylim([0 3]);
title(char(strcat('Normalized diffusion map prototype spectra, K=',num2str(k),', colored by log age')));


% plot normalized prototype spectra, colored by log Z
CFcolor = ceil((log10(Z(idxCF))-min(log10(protometals)))/(max(log10(protometals))-min(log10(protometals)))*63)+1;
dmapcolor = ceil((log10(protometals)-min(log10(protometals)))/(max(log10(protometals))-min(log10(protometals)))*63)+1;

figure, % prototype spectra, colored by log Z
subplot(2,1,1)
for ii=1:45,
  plot(lam,ssps(idxCF(ii),:),'Color',cmap(CFcolor(ii),:),'LineWidth',.01)
  hold on
end
xlabel('lambda'); ylabel('Normalized flux');ylim([0 3]);
title('Normalized Cid Fernandes et al.(2005) prototype spectra, K=45, colored by log Z');
subplot(2,1,2)
for ii=1:k,
  plot(lam,protospecnorm(ii,:),'Color',cmap(dmapcolor(ii),:),'LineWidth',.01)
  hold on
end
xlabel('lambda'); ylabel('Normalized flux'); ylim([0 3]);
title(char(strcat('Normalized diffusion map prototype spectra, K=',num2str(k),', colored by log Z')));


