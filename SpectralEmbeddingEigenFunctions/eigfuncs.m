function [dd2,uu2] = eigfuncs(xx,NUM_EVECS,figs) 
% adapted from fergus' demo: more pictures, more general, binning slightly
% more stable (Wsquig), also sigma not constant, except for exact solution
 
nPoints = prod(size(xx));
ndata = size(xx);

if nargin < 2
    NUM_EVECS = 15;
end
if nargin < 3
    figs = 0; % don't plot figures!
end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Classify with Eigenfunctions!
 
 % compute approximate eigenvectors using eigenfunction approach
 EPSILON = 0% .01;
 %NUM_EVECS = 15;
 SIGMA = .7;
 [dd2,uu2] = eigenfunctions_image(xx,SIGMA,NUM_EVECS,EPSILON,figs);
 dd2d = diag(dd2);

 % paint generalized eigenfunctions onto data... 
h4 = figure; 
nplots = floor(sqrt(NUM_EVECS) + 1);
for k = 1:NUM_EVECS
    subplot(nplots,nplots,k);
    imagesc(reshape(uu2(:,k),ndata));
    title(['eigfunc_',num2str(k),' ,sigma_',num2str(k),' = ',num2str(dd2d(k))],'FontSize',15);
end
 
 

 