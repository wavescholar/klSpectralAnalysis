clear all;
set(0,'DefaultFigureWindowStyle','docked')
%set(0,'DefaultFigureWindowStyle','normal')
global figPath
figPath ='D:\klSpectralAnalysis\output\';

%Write Generate Data Function
mu1 = [-20 0];
sigma1 = [3 .4; .4 3];
mu2 = [-15 0];
sigma2 = [3 0; 0 3];
mu3 = [-4 0];
sigma3 = [3 0; 0 1.5];
mu4 = [4 0];
sigma4 = [2 0; 0 1];

X1=[mvnrnd(mu1,sigma1,3200)];
X2=[mvnrnd(mu2,sigma2,3200)];
X3=[mvnrnd(mu3,sigma3,3200)];
X4 =[ mvnrnd(mu4,sigma4,3200)];
X=[X1;X2;X3;X4];
figure ;
scatter(X1(:,1),X1(:,2),'b','.');hold on;
scatter(X2(:,1),X2(:,2),'g','.');
scatter(X3(:,1),X3(:,2),'y','.');
scatter(X4(:,1),X4(:,2),'r','.');
data =X;

%We have 2 D data. Lets make some affinity matrices.
nrand = 1000;
minneighbors = 5;
nframes = size(data,1)
rr = floor(nframes*rand(nrand,1)) + 1;
mindists = zeros(nrand,1);
for k = 1:nrand
    kk = rr(k);
    %    qdist = sum((repmat(qaa(kk,:),nframes,1)-qaa(1:nframes,:)).^2,2);
    qdist = sum((ones(nframes,1)*data(kk,:)-data).^2,2);
    qdists = sort(qdist);
    mindists(k) = qdists(minneighbors);
end
maxrad = median(mindists)
[Aloc ii jj maxrad] = makeAffinity_local_cerno(data,minneighbors,maxrad);

%For testing in KL C++ MAtrix Library. 
%spy(Aloc);
L2=laplacian(Aloc);
csvwrite(['L_',num2str(size(X,1)),'.txt'],full(L2));
[V,D] = eigs(L2,16,'SA');
for j=1:16
    csvwrite(['EV_',num2str(j-1),'_L_',num2str(size(X,1)),'.txt'],V(:,j));
end

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
h =   figure;%('Visible','off'); 
set(h,'Position',[ 0 0 1724 945]);
for k = 1:nhh
    subplot(ns,ns,k);
    spy(hh{k}.A,'k'); formatfig;
    title(['level ',num2str(k),', ',num2str(size(hh{k}.A,1)),' clusters'],'FontSize',16);
    set(gca,'FontSize',12);
end
%export_fig(h,[ figPath,'AllTiles_hh_A.jpg']);
print(gcf,'-dpdf','-r600',[ figPath,'AllTiles_hh_A.pdf']);
close all;

CM=colormap(lines);
for level=7:9
    figure;%h=figure('Visible','off'); 
    hold on;
    %set(h,'Position',[ 32 32 1820 980]);
    for i=1:size(data,1)
        clusterNumber =hh{level}.cl(i);
        scatter(data(i,1),data(i,2),'Marker','.','MarkerEdgeColor' ,CM(clusterNumber,:));
    end
    print(gcf,'-dpdf','-r600',[ figPath,'class',num2str(level),'.pdf']);
    %export_fig(h,[ figPath,'class',num2str(level),'.jpg']);
    
    close all;
end
















clear all;
N=3*   2594;
x = (1:N);%randperm(N);
gs = N/3;
G1 = x(1:gs);
G2 = x(gs+1:2*gs);
G3= x(2*gs+1:3*gs);
p_G1 = 0.000009;
p_G2 = 0.000004;
p_G3 = 0.000005;
p_Inbetween = 0.0000009;
A(G1, G1) = rand(gs,gs) < p_G1;
A(G2, G2) = rand(gs,gs) < p_G2;
A(G3, G3) = rand(gs,gs) < p_G3;
[n1,m1]=size( A(G1, G2));
[n2,m2]=size( A(G2, G3));
[n3,m3]=size( A(G3, G1));
B_1=rand(gs, gs);B_2=rand(gs, gs);B_3=rand(gs, gs);
 A(G1, G2) = B_1 < p_Inbetween;
 A(G2, G3) = B_2 < p_Inbetween;
 A(G1, G3) = B_3 < p_Inbetween;
%clear B_1 B_2 B_3;
A = triu(A,1);
A = A + A';
spy(A);
subsp = 2;
nrounds = 50;
hh = svdms_msm(A,subsp,nrounds);



%L = del2(A);
L2=laplacian(Aloc);%
Ic = adjacency2incidence(A);
%L3 = Ic*Ic';%compute_combinatorial_laplacian(A);
[V D] = eigs(L2, 3, 'SA');sigma=eig(A);hist(sigma,50);
D(2,2);plot(V(:,2), '.-');
plot(sort(V(:,2)), '.-');
plot(V(:,3), '.-');
plot(sort(V(:,3)), '.-');
[ignore p] = sort(V(:,2));
spy(A(p,p));

XY_ga = gursoy_atun_layout(sparse(A*1.0),'topology','circle');
XY_kk = kamada_kawai_spring_layout(sparse(A*1.0));
XY_fr=fruchterman_reingold_force_directed_layout(sparse(A*1.0),'initial_temp',50);
wgPlot(A,XY_fr); %This is the best format for cut display
gplot(A,XY_fr);