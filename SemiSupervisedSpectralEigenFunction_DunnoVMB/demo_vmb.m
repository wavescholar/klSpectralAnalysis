 
% adapted from fergus' demo: more pictures, more general, binning slightly
% more stable (Wsquig), also sigma not constant, except for exact solution
 
 clear; close all;
 
 %%% user parameters
 NUM_EVECS = 5; % how many eigenvectors/eigenfunctions to use
 SIGMA = 0.2; % controls affinity in graph Laplacian

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% generate data
 
 square = 0;
 if square
     %%% generate some random data.
     xrange=[0:0.045:1];
     yrange=[0:0.045:2];
     [xx,yy]=meshgrid(xrange,yrange);
     xx=[xx(:) yy(:)];
     xx=xx+0.01*randn(size(xx));
     nPoints = size(xx,1);
     norig = nPoints;
 else
     % TWO GAUSSIANS:
     xrange=[0:0.045:1];
     yrange=[0:0.045:2];
     [xx,yy]=meshgrid(xrange,yrange);
     xx=[xx(:) yy(:)];
     nPoints = size(xx,1)
     halfn = floor(nPoints/2);
     mu1 = [2 1];
     Sigma1 = [1 .5; .5 2];%[1 0; 0 25]; %[1 .5; .5 2];
     R1 = chol(Sigma1);
     xx(1:halfn,:) = repmat(mu1,halfn,1) + randn(halfn,2)*R1;
     mu2 = [10 1];
     Sigma2 =  [1 .5; .5 2];%[1 0; 0 25]; %[1 .1; 1 .95]; %[1 .5; .5 2];
     R2 = chol(Sigma2);
     xx(halfn+1:end,:) = repmat(mu2,halfn+1,1) + randn(halfn+1,2)*R2;
     norig = nPoints;
     % Instead of regularizing histogram, add a uniform distribution of points
     % to data:
     min(xx)
     max(xx)
     % add points to data now, so that there is data everywhere (instead of
     % regularizing histogram)
     nreg = round(nPoints/10);
     basisxx(:,1) = min(xx(:,1)) + (max(xx(:,1)) - min(xx(:,1)))*rand(nreg,1);
     basisxx(:,2) = min(xx(:,2)) + (max(xx(:,2)) - min(xx(:,2)))*rand(nreg,1);
     data = [xx; basisxx];
     xx = data;
     nPoints  = size(xx,1)
 end
 
figure; subplot(2,2,1); scatter(xx(:,1),xx(:,2)); 
xlabel('x1','FontSize',20);
ylabel('x2','FontSize',20);
subplot(2,2,3); [n s] = hist(xx(:,1),max(nPoints/20,4)); bar(s,n/sum(n)); title('x1','FontSize',20);
subplot(2,2,2); [n s] = hist(xx(:,2),max(nPoints/20,4)); barh(s,n/sum(n)); title('x2','FontSize',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% option: do PCA
[co pc] = princomp(xx); 
figure; plot(pc(:,1),pc(:,2),'rx')
xx = pc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% select points to label!
nDims = size(xx,2);
minxx = zeros(nDims,1);
maxx = zeros(nDims,1);
cutperc = .05; % go with .05, bc 95% of data is within 3sigma of center for gaussian distributed...
for k = 1:nDims
    sortx = sort(xx(:,k));
    minxx(k) = sortx(round(cutperc*nPoints));
    maxx(k) = sortx(round((1-cutperc)*nPoints));
end

nlabels = 2;
is = zeros(nlabels,1);
[dummy,is(1)]=min((xx(:,1)-minxx(1)).^2+(xx(:,2)-minxx(2)).^2);
[dummy,is(2)]=min((xx(:,1)-maxx(1)).^2+(xx(:,2)-maxx(2)).^2);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Plot data and labeled points
 h=figure;
 subplot(1,4,1);
 hold off;
 plot(xx(:,1),xx(:,2),'x');
 hold on;
 plot(xx(is(1),1),xx(is(1),2),'ro','Markersize',8);
 plot(xx(is(2),1),xx(is(2),2),'go','Markersize',8);
 axis equal;
 axis off;
 title('Data','Fontsize',16);
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% setup weights on datapoints (lambda) and labels (y)
lambda=zeros(nPoints,1);
y = zeros(nPoints,1);
y(is(1))=+1;
lambda(is(1))=1000;
y(is(2))=-1;
lambda(is(2))=1000;
 % build diagonal Lambda matrix
 Lambda=diag(lambda);
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Compute normalized graph Laplacian
 dists = dist_mat(xx,xx);
 SIGMA = round(median(sqrt(dists(:))))
 W=exp(-0.5*dists.^2/SIGMA^2);
 D=diag(sum(W));
 L=D-W;
 L=inv(sqrt(D))*L*inv(sqrt(D)); % graph laplacian made with data-dependent sigma


 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Classify with Eigenfunctions!
 
 % compute approximate eigenvectors using eigenfunction approach
 EPSILON = 0% .01;
 [dd2,uu2] = eigenfunctions_vmb(xx,SIGMA,NUM_EVECS,EPSILON);
 dd2d = diag(dd2);

 % paint generalized eigenfunctions onto data... 
h4 = figure; 
nplots = min(NUM_EVECS,5);
for k = 1:nplots
    subplot(nplots,2,2*k);
    paintArvind(uu2(:,k),xx);
    title(['eigfunc_',num2str(k),' ,sigma_',num2str(k),' = ',num2str(dd2d(k))],'FontSize',15);
end
 
 
 %%% now solve for coefficients in NUM_EVECS x NUM_EVECS linear system
 %%% this is eqn. 1 in the NIPS paper (but with approx. eigenvectors/values).
 alpha2=(dd2 +uu2'*Lambda*uu2)\(uu2'*Lambda*y);
 f_efunc=uu2*alpha2;
 
 % plot out
 figure(h);
 subplot(1,4,4);
 Ipos=find(f_efunc>0);
 Ineg=find(f_efunc<0);
 hold off;
 plot(xx(Ipos,1),xx(Ipos,2),'ro');
 hold on;
 plot(xx(Ineg,1),xx(Ineg,2),'go');
 axis equal;
 axis off;
 title('Eigenfunction','Fontsize',16);


 %%%%%%%%%%%%%%%%%%%%%%%%%%%5
 %% Compare eigenfunction solution to 
 % 1) solution with least-squares problem
 % 2) solution with eigenvectors. 
 
 
 % use a small sigma for the exact solution (just seems to work better,
 % and Fergus uses SIGMA = 0.2)
%  SIGMA = .2% round(median(sqrt(dists(:))))
%  W=exp(-0.5*dists.^2/SIGMA^2);
%  D=diag(sum(W));
%  L=D-W;
%  L=inv(sqrt(D))*L*inv(sqrt(D));
 L_lsq = L; % graph laplacian made with SIGMA = .2
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Solve full least-squares system for labels
 f_lsq=(L_lsq+Lambda)\(Lambda*y);
 
 %% plot out
 figure(h);
 subplot(1,4,2);
 Ipos=find(f_lsq>0);
 Ineg=find(f_lsq<0);
 hold off;
 plot(xx(Ipos,1),xx(Ipos,2),'ro');
 hold on;
 plot(xx(Ineg,1),xx(Ineg,2),'go');
 axis equal;
 axis off;
 title('Exact LSQ','Fontsize',16);
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Find eigenvectors uu of Laplacian and solve 
 %% within this basis: f = uu * alpha 
 
 %%% compute exact eigenvectors/values of graph Laplacian
 [uu,ss,vv]=svd(L);
 %uu=fliplr(uu(:,end-NUM_EVECS:end-1));
 uu=fliplr(uu(:,end-NUM_EVECS:end));
 dd = diag(ss);
 % dd = dd(end-NUM_EVECS:end-1); 
 dd = dd(end-NUM_EVECS:end); 
 dd =flipud(dd);
 uu=inv(sqrt(D))*uu;

 % vmb: paint generalized eigenvectors onto data...
figure(h4); 
nplots = min(NUM_EVECS,5);
for k = 1:nplots
    subplot(nplots,2,2*k-1);
    paintArvind(uu(:,k),xx);
    title(['eigvec_',num2str(k),' ,sigma_',num2str(k),' = ',num2str(dd(k))],'FontSize',15);
end



 %%% now solve for coefficients in NUM_EVECS x NUM_EVECS linear system 
 %%% this is eqn. 1 in the NIPS paper.
 alpha=(uu'*L*uu+uu'*Lambda*uu)\(uu'*Lambda*y);
 f_evec=uu*alpha;
 
 % plot out 
 figure(h);
 subplot(1,4,3);
 Ipos=find(f_evec>0);
 Ineg=find(f_evec<0);
 hold off;
 plot(xx(Ipos,1),xx(Ipos,2),'ro');
 hold on;
 plot(xx(Ineg,1),xx(Ineg,2),'go');
 axis equal;
 axis off;
 title('Eigenvector','Fontsize',16);
 
 
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% Now compare exact/approximate eigenvectors
 h2=figure;
 nplots = min(NUM_EVECS,5);
 for a=1:nplots%NUM_EVECS
   % exact eigenvectors
   subplot(NUM_EVECS,2,a*2-1);
   
   % put eigvecs and eigfuncs on same scale!
   eigvecs = reshape(uu(1:norig,a),[length(yrange) length(xrange)]);
   eigfuncs = reshape(uu2(1:norig,a),[length(yrange) length(xrange)]);
   m = min([min(eigvecs(:)) min(eigfuncs(:))]);
   M = max([max(eigvecs(:)) max(eigfuncs(:))]);
   
   imagesc(eigvecs,[m M]); colorbar;
   title(['eigval = ',num2str(dd(a))],'FontSize',16);
   if (a==1)
     title(['Exact eigenVECTORs, eigval = ',num2str(dd(a))],'FontSize',16);
   end
   axis equal;  axis off;
   
  
   % approximate eigenvectors
   subplot(NUM_EVECS,2,a*2);
   imagesc(eigfuncs,[m M]);colorbar;
   title(['eigval = ',num2str(dd2(a,a))],'FontSize',16);
   if (a==1)
     title(['approx. eigenFUNCTIONs, eigval = ',num2str(dd2(a,a))],'FontSize',16);
   end
   axis equal; axis off
   
 end
 set(h2,'Position',[ 127   753   861   717]);
 
figure;
for k = 1:min(NUM_EVECS,9);
    subplot(3,3,k);
    plot(uu(:,k)); hold on; plot(uu2(:,k),'r')
    title(['eigvec_',num2str(k),' blue, eigfunc_',num2str(k),' red'],'FontSize',14);
end


