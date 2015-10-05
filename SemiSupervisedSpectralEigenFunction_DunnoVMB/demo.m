 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Demo for eigenfunction method for computing approximate 
 % eigenvectors of graph Laplacian. 
 %
 % This function runs three different variant of semi-supervised
 % learning using (normalized) graph Laplacians:
 % 1. Exact least-squares solution.
 % 2. Solution using smallest few exact eigenvectors of graph Laplacian.
 % 3. Solution using smallest few approximate eigenvectors of graph
 % Laplacian. The approximate eigenvectors are computed from the
 % numerical eigenfunctions method described in the NIPS 2009 paper:
 % "Semi-supervised learning in gigantic image collections" 
 % by R. Fergus, Y. Weiss and A. Torralba. 
 %
 % The code generates some toy data and a pair of labels, runs each of
 % the 3 methods above on it and then visualizes the exact and
 % approximate eigenvectors/values. 
 %
 % The exact least-squares is a little bit unstable, so you'll see the
 % decision boundary jumping around if you run the demo a few times. 
 %
 % Note that methods 1 and 2 are polynomial in the number of datapoints,
 % while method 3 is linear in the number datapoints (but approximate).
 %
 % The implementation of the numerical eigenfunctions approach is general
 % and can be applied to real data (it is a cleaned up version of the
 % code used to generate the results in the NIPS paper). Note that it
 % does assume seperability of the data, so you might want to preprocess
 % it first with PCA (something that this demo does NOT do). 
 %
 % If you do use this code in your paper, please cite our paper. For your
 % convenience, here is the bibtex: 
 %
 % @inproceedings{FergusWeissTorralbaNIPS09,
 %  author = "Fergus, R. and Weiss, Y. and Torralba, A.",
 %  title = "Semi-supervised Learning in Gigantic Image Collections",
 %  booktitle = "NIPS",
 %  year = "2009"
 % } 
 %
 % Version 1.0. Rob Fergus & Yair Weiss. 11/30/09.  
  
 
 clear; close all;
 
 %%% user parameters
 NUM_EVECS = 5; % how many eigenvectors/eigenfunctions to use
 SIGMA = 0.2; % controls affinity in graph Laplacian

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Setup data, labels and build graph Laplacian
 
 %%% generate some random data. 
 xrange=[0:0.05:1];
 yrange=[0:0.05:2];
 [xx,yy]=meshgrid(xrange,yrange);
 xx=[xx(:) yy(:)];
 xx=xx+0.01*randn(size(xx));
 nPoints = size(xx,1);  

 %%% Compute normalized graph Laplacian
 W=exp(-0.5*dist_mat(xx,xx).^2/SIGMA^2);
 D=diag(sum(W));
 L=D-W;
 L=inv(sqrt(D))*L*inv(sqrt(D));

 %%% setup weights on datapoints (lambda) and labels (y)
 lambda=zeros(nPoints,1);
 y = zeros(nPoints,1);

 % one positive example at (0,0)
 % and one negative example at (1,2)
 [dummy,i1]=min(xx(:,1).^2+xx(:,2).^2);
 y(i1)=+1;
 lambda(i1)=1000;
 [dummy,i2]=min((xx(:,1)-1).^2+(xx(:,2)-2).^2);
 y(i2)=-1;
 lambda(i2)=1000;

 % build diagonal Lambda matrix
 Lambda=diag(lambda);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Plot data and labeled points
 h=figure;
 subplot(1,4,1);
 hold off;
 plot(xx(:,1),xx(:,2),'x');
 hold on;
 plot(xx(i1,1),xx(i1,2),'ro','Markersize',8);
 plot(xx(i2,1),xx(i2,2),'go','Markersize',8);
 axis equal;
 axis off;
 title('Data','Fontsize',16);
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Solve full least-squares system for labels
 f_lsq=(L+Lambda)\(Lambda*y);
 
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
 uu=fliplr(uu(:,end-NUM_EVECS:end-1));
 dd = diag(ss);
 dd = dd(end-NUM_EVECS:end-1); dd =flipud(dd);
 uu=inv(sqrt(D))*uu;

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
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Now try approximate eigenvectors, as found
 %% by eigenfunction approach
 
 % compute approximate eigenvectors using eigenfunction approach
 [dd2,uu2] = eigenfunctions(xx,SIGMA,NUM_EVECS);

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


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% Now compare exact/approximate eigenvectors
 h2=figure;
 for a=1:NUM_EVECS
   % exact eigenvectors
   subplot(NUM_EVECS,3,a*3-2);
   imagesc(reshape(uu(:,a),[length(yrange) length(xrange)]));
   if (a==1)
     title(sprintf('Exact\n Eigenvectors'),'Fontsize',16);
   end
   axis equal;  axis off;
   
   % eigenvalues
   subplot(NUM_EVECS,3,a*3-1);
   text(0.1,0.5,[num2str(dd(a),'%0.4f'),' -- ',num2str(dd2(a,a),'%0.4f')],'Fontsize',16);
   axis off
   if (a==1)
     title(sprintf('Exact (L) & \nApproximate (R)\n Eigenvalues'),'Fontsize',16);
   end
   
   % approximate eigenvectors
   subplot(NUM_EVECS,3,a*3);
   imagesc(reshape(uu2(:,a),[length(yrange) length(xrange)]));
   if (a==1)
     title(sprintf('Approximate\n Eigenvectors'),'Fontsize',16);
   end
   axis equal; axis off
   
  
 end
 set(h2,'Position',[ 127   753   861   717]);
 
