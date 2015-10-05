clear all;
mu1 = [-20 0];
sigma1 = [3 .4; .4 3];

mu2 = [-15 0];
sigma2 = [3 0; 0 3];

mu3 = [-4 0];
sigma3 = [3 0; 0 1.5];

mu4 = [4 0];
sigma4 = [2 0; 0 1];

X1=[mvnrnd(mu1,sigma1,1000)]
X2=[mvnrnd(mu2,sigma2,1000)];
X3=[mvnrnd(mu3,sigma3,1000)];
X4=[ mvnrnd(mu4,sigma4,1000)];
X=[X1;X2;X3;X4];

figure ;hold on;
scatter(X1(:,1),X1(:,2),'b','.');
scatter(X2(:,1),X2(:,2),'g','.');
scatter(X3(:,1),X3(:,2),'y','.');
scatter(X4(:,1),X4(:,2),'r','.');
fileName= 'Data.pdf';
print(gcf,'-dpdf',fileName);close gcf; 

%%
%PCA SVD 
[coeff,score,latent] = pca(X1);
[Ve,De] =eigs(cov(X1))
[U,S,V] = svd(X1);

figure ;hold on;
gm = gmdistribution.fit(X1,1);
ezcontour(@(x,y)pdf(gm,[x y]),[-30 10],[-10 10],160);
covX1 = cov(X1);
quiver(mu1(1),mu1(2) ,covX1(1,1),covX1(2,1),'b');
quiver(mu1(1),mu1(2) ,covX1(1,2),covX1(2,2),'b');

gm = gmdistribution.fit(X2,1);
ezcontour(@(x,y)pdf(gm,[x y]),[-30 10],[-10 10],160);
covX2 = cov(X2);
quiver(mu2(1),mu2(2) ,covX2(1,1),covX2(2,1),'g');
quiver(mu2(1),mu2(2) ,covX2(1,2),covX2(2,2),'g');

gm = gmdistribution.fit(X3,1);
ezcontour(@(x,y)pdf(gm,[x y]),[-30 10],[-10 10],160);
covX3 = cov(X3);
quiver(mu3(1),mu3(2) ,covX3(1,1),covX3(2,1),'y');
quiver(mu3(1),mu3(2) ,covX3(1,2),covX3(2,2),'y');

gm = gmdistribution.fit(X4,1);
ezcontour(@(x,y)pdf(gm,[x y]),[-30 10],[-10 10],630);
covX4 = cov(X4);
quiver(mu4(1),mu4(2) ,covX4(1,1),covX4(2,1),'r');
quiver(mu4(1),mu4(2) ,covX4(1,2),covX4(2,2),'r');
fileName= 'Model.pdf';
print(gcf,'-dpdf',fileName);close all;

%%
% [U,S,V] = svd(X1);
% [Ve,De] =eigs(sigma1);
% mean(X1)
% std(X1)
% 
% %Standardize - this is a faster way to do it for large arrays
% C = bsxfun(@minus,X1,mean(X1));
% Xstd = bsxfun(@rdivide,C,std(X1));
% mean(Xstd)
% std(Xstd)
% [coeff_std,score_std,latent_std] = pca(Xstd);
% [Ustd,Sstd,Vstd] = svd(Xstd);
% 
% figure ;hold on;
% scatter(X1(:,1),X1(:,2),'b','.');
% quiver(mu1(1),mu1(2) ,V(1,1),V(2,1),'r');
% quiver(mu1(1),mu1(2) ,latent(1,1),latent(1,2),'r');
% 
% 
% figure ;hold on;gm = gmdistribution.fit(X1,1);
% ezcontour(@(x,y)pdf(gm,[x y]),[-30 0],[-5 5],60);
% 
% scatter(X2(:,1),X2(:,2),'g','.');
% scatter(X3(:,1),X3(:,2),'y','.');
% scatter(X4(:,1),X4(:,2),'r','.');
% fileName= 'Data.pdf';
% print(gcf,'-dpdf',fileName);close gcf; 

% hold on;gm = gmdistribution.fit(X,4);
% ezcontour(@(x,y)pdf(gm,[x y]),[-30 10],[-10 10],60);
% fileName= 'DataFitToMixtureOf4Gaussians.pdf';
% print(gcf,'-dpdf',fileName);close gcf; 


%%
nrand = 1000;
minneighbors = 5;
nframes = size(X,1)
rr = floor(nframes*rand(nrand,1)) + 1;
mindists = zeros(nrand,1);
for k = 1:nrand
    kk = rr(k);
    qdist = sum((ones(nframes,1)*X(kk,:)-X).^2,2);
    qdists = sort(qdist);
    mindists(k) = qdists(minneighbors);
end

maxrad = median(mindists)
[Aloc ii jj maxrad] = makeAffinity_local_ZelnikManorPerona(X,minneighbors,maxrad);

figure ; spy(Aloc);
fileName= 'spy_A_Zelnik_Manor.pdf';
print(gcf,'-dpdf',fileName);close gcf; 

L2=laplacian(Aloc);
[V,D] = eigs(L2,16,'SA');
figure ; plot(sort(V(:,1)), 'c');hold on;
plot(sort(V(:,2)), 'b');
plot(sort(V(:,3)), 'r');
plot(sort(V(:,4)), 'g');
title('First Four Eigenvectors of Graph Laplacian')
legend('EV 0','EV1','EV 2','EV 3');
print(gcf,'-dpdf','-r600',[ 'EigenVecsOfGraphLaplacian.pdf']);
close gcf;

subsp = 2;
nrounds = 50;
hh = svdms_msm(Aloc,subsp,nrounds);
if size(hh{end}.A,1) == 1
    hh = hh(1:end-1);
end
nhh = size(hh,2);
if floor(sqrt(nhh)) == sqrt(nhh)
    ns = sqrt(nhh)
else
    ns = floor(sqrt(nhh)) + 1
end
h =   figure;
set(h,'Position',[ 0 0 1724 945]);
for k = 1:nhh
    subplot(ns,ns,k);
    spy(hh{k}.A,'k'); formatfig;
    title(['level ',num2str(k),', ',num2str(size(hh{k}.A,1)),' clusters'],'FontSize',16);
    set(gca,'FontSize',12);
end
print(gcf,'-dpdf','-r600',[ 'HeirarchicalEigensolver_AffinityMatrics.pdf']);
close all;

CM=colormap(lines);
for level=nhh-4:nhh
    figure;
    hold on;
    for i=1:size(X,1)
        clusterNumber =hh{level}.cl(i);
        scatter(X(i,1),X(i,2),'Marker','.','MarkerEdgeColor' ,CM(clusterNumber,:));
    end
    title(['EigenCut Class Labels Level ',num2str(level)]); 
    print(gcf,'-dpdf','-r600',[ 'EigenCutClassLabels_Level',num2str(level),'.pdf']);
  close all;
end

%%Diffusion Maps

riskClass=zeros(1,size(X,1));
riskClass(1:1000)=0;
riskClass(1001:2000)=1;
riskClass(2001:3000)=2;
riskClass(3001:4000)=3;


%Standardize - this is a faster way to do it for large arrays
C = bsxfun(@minus,X,mean(X));
X = bsxfun(@rdivide,C,std(X));
mean(X)
std(X)

graphDiffusion_opts = struct( ...
    'kNN'           , 100, ...                          % How many nearest neighbors each point is connected to
    'kNNAutotune'   , 15, ...                           % Local scale parameter = distance to the kNNAutotune-th nearest neighbor
    'Normalization' , 'smarkov', ...                    % How to normalize diffusion matrix
    'kEigenVecs'    , 500 );                            % How many eigenvectors to compute
G = GraphDiffusion(X', 0, graphDiffusion_opts);

figure;plot(G.EigenVals);title('Plot of diffusion eigenvalues');
fileName= 'DiffusionMapEigenValues.pdf';
print(gcf,'-dpdf',fileName);close gcf; 

% Risk as a function on diffusion embedding
figure;
for k = 2:6,
    figure ;
    gscatter(G.EigenVecs(:,k),G.EigenVecs(:,k+1),riskClass,[],[],10);
    legend('Class 1','Class 2','Class3','Class 4');
    title(sprintf('Coordinates %d,%d, color=risk ',k,k+1));
    fileName= [sprintf('DiffusionEigenVectors%d%d',k,k+1),'.pdf'];
print(gcf,'-dpdf',fileName);close gcf; 
close gcf;
end;

% One can learn to discriminate between RISK using eigenvectors.
% Approximate a disciminating function between Risk with diffusion eigenvectors
F_coeffs = riskClass*G.EigenVecs;
F=F_coeffs*G.EigenVecs';                               % Expand function on eigenvectors

figure;plot(F_coeffs);title('Coefficient of regression function');
fileName= 'CoefficientOfRegressionFunction.pdf';
print(gcf,'-dpdf',fileName);close gcf; 

figure;plot(cumsum(F_coeffs.^2)/norm(F));
title('Energy capture as a function of number of coefficients kept');
fileName= 'EnergyCaptureAsFunctionOfNumberOfCoefficientsKept.pdf'
print(gcf,'-dpdf',fileName);close gcf; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%Now Embed In 3D
clear all;
mu1 = [-20 0];
sigma1 = [3 .4; .4 3];

mu2 = [-15 0];
sigma2 = [3 0; 0 3];

mu3 = [-4 0];
sigma3 = [3 0; 0 1.5];

mu4 = [4 0];
sigma4 = [2 0; 0 1];

X1=[mvnrnd(mu1,sigma1,1000)]
X2=[mvnrnd(mu2,sigma2,1000)];
X3=[mvnrnd(mu3,sigma3,1000)];
X4=[ mvnrnd(mu4,sigma4,1000)];
X=[X1;X2;X3;X4];
X3 = [X3,atan(X3(:,1))];
X2 = [X2,atan(X2(:,1))];
X1 = [X1,atan(X1(:,1))];
X4 = [X4,atan(X4(:,1))];
X=[X1;X2;X3;X4];

figure ;
scatter3(X1(:,1),X1(:,2),X1(:,3),'b','.');hold on;
scatter3(X2(:,1),X2(:,2),X2(:,3),'g','.');
scatter3(X3(:,1),X3(:,2),X3(:,3),'y','.');
scatter3(X4(:,1),X4(:,2),X4(:,3),'r','.');
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
[x,y] = meshgrid(-30:.5:10,-10:.5:10);
z = atan(x);
h=surf(x,y,z,'EdgeColor',[.7,.7,.7]);alpha('color')
colormap(bone) 
fileName= '3DData.pdf';
print(gcf,'-dpdf',fileName);close gcf; 
% 
% hold on;
% gm = gmdistribution.fit(X,4);
% ezcontour(@(x,y,z)pdf(gm,[x y,z]),[-30 10],[-10 10],[-1,1],60);
% fileName= '3DDataFitToMixtureOf4Gaussians.pdf';
% print(gcf,'-dpdf',fileName);close gcf; 

nrand = 1000;
minneighbors = 5;
nframes = size(X,1)
rr = floor(nframes*rand(nrand,1)) + 1;
mindists = zeros(nrand,1);
for k = 1:nrand
    kk = rr(k);
    qdist = sum((ones(nframes,1)*X(kk,:)-X).^2,2);
    qdists = sort(qdist);
    mindists(k) = qdists(minneighbors);
end

maxrad = median(mindists)
[Aloc ii jj maxrad] = makeAffinity_local_ZelnikManorPerona(X,minneighbors,maxrad);

figure ; spy(Aloc);
fileName= '3Data_Spy_A_Zelnik_Manor.pdf';
print(gcf,'-dpdf',fileName);close gcf; 

L2=laplacian(Aloc);
[V,D] = eigs(L2,16,'SA');
figure ; plot(sort(V(:,1)), 'c');hold on;
plot(sort(V(:,2)), 'b');
plot(sort(V(:,3)), 'r');
plot(sort(V(:,4)), 'g');
title('First Four Eigenvectors of Graph Laplacian')
legend('EV 0','EV1','EV 2','EV 3');
print(gcf,'-dpdf','-r600',[ '3DData_EigenVecsOfGraphLaplacian.pdf']);
close gcf;

subsp = 2;
nrounds = 50;
hh = svdms_msm(Aloc,subsp,nrounds);
if size(hh{end}.A,1) == 1
    hh = hh(1:end-1);
end
nhh = size(hh,2);
if floor(sqrt(nhh)) == sqrt(nhh)
    ns = sqrt(nhh)
else
    ns = floor(sqrt(nhh)) + 1
end
h =   figure;
set(h,'Position',[ 0 0 1724 945]);
for k = 1:nhh
    subplot(ns,ns,k);
    spy(hh{k}.A,'k'); formatfig;
    title(['level ',num2str(k),', ',num2str(size(hh{k}.A,1)),' clusters'],'FontSize',16);
    set(gca,'FontSize',12);
end
print(gcf,'-dpdf','-r600',[ '3DData_HeirarchicalEigensolver_AffinityMatrics.pdf']);
close all;

CM=colormap(lines);
for level=nhh-2:nhh
    figure;
    hold on;
    for i=1:size(X,1)
        clusterNumber =hh{level}.cl(i);
        scatter3(X(i,1),X(i,2),X(i,3),'Marker','.','MarkerEdgeColor' ,CM(clusterNumber,:));
    end
    title(['EigenCut Class Labels Level ',num2str(level)]); 
    print(gcf,'-dpdf','-r600',[ '3DData_EigenCutClassLabels_Level',num2str(level),'.pdf']);
  close all;
end

%Diffusion Maps

riskClass=zeros(1,size(X,1));
riskClass(1:1000)=0;
riskClass(1001:2000)=1;
riskClass(2001:3000)=2;
riskClass(3001:4000)=3;


%Standardize - this is a faster way to do it for large arrays
C = bsxfun(@minus,X,mean(X));
X = bsxfun(@rdivide,C,std(X));
mean(X)
std(X)

graphDiffusion_opts = struct( ...
    'kNN'           , 100, ...                          % How many nearest neighbors each point is connected to
    'kNNAutotune'   , 15, ...                           % Local scale parameter = distance to the kNNAutotune-th nearest neighbor
    'Normalization' , 'smarkov', ...                    % How to normalize diffusion matrix
    'kEigenVecs'    , 500 );                            % How many eigenvectors to compute
G = GraphDiffusion(X', 0, graphDiffusion_opts);

figure;plot(G.EigenVals);title('Plot of diffusion eigenvalues');
fileName= '3DData_DiffusionMapEigenValues.pdf';
print(gcf,'-dpdf',fileName);close gcf; 

% Risk as a function on diffusion embedding
figure;
for k = 2:6,
    figure ;
    gscatter(G.EigenVecs(:,k),G.EigenVecs(:,k+1),riskClass,[],[],10);
    legend('Class 1','Class 2','Class3','Class 4');
    title(sprintf('Coordinates %d,%d, color=risk ',k,k+1));
    fileName= [sprintf('3DData_DiffusionEigenVectors%d%d',k,k+1),'.pdf'];
print(gcf,'-dpdf',fileName);close gcf; 
close gcf;
end;

% One can learn to discriminate between RISK using eigenvectors.
% Approximate a disciminating function between Risk with diffusion eigenvectors
F_coeffs = riskClass*G.EigenVecs;
F=F_coeffs*G.EigenVecs';                               % Expand function on eigenvectors

figure;plot(F_coeffs);title('Coefficient of regression function');
fileName= '3DData_CoefficientOfRegressionFunction.pdf';
print(gcf,'-dpdf',fileName);close gcf; 

figure;plot(cumsum(F_coeffs.^2)/norm(F));
title('Energy capture as a function of number of coefficients kept');
fileName= 'EnergyCaptureAsFunctionOfNumberOfCoefficientsKept.pdf'
print(gcf,'-dpdf',fileName);close gcf; 















