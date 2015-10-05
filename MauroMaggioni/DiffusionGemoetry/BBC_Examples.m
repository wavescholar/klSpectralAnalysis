X=[cos(linspace(0,2*pi-2/500*pi,500));sin(linspace(0,2*pi-2/500*pi,500))];

figure ;plot(X(1,:),X(2,:))
graphDiffusion_opts = struct( ...
    'kNN'           , 100, ...                          % How many nearest neighbors each point is connected to
    'kNNAutotune'   , 15, ...                           % Local scale parameter = distance to the kNNAutotune-th nearest neighbor
    'Normalization' , 'smarkov', ...                    % How to normalize diffusion matrix
    'kEigenVecs'    , 500 );                            % How many eigenvectors to compute
G = GraphDiffusion(allTileData_norm, 0, graphDiffusion_opts);

graphDiffusion_opts = struct( ...
    'kNN'           , 100, ...                          % How many nearest neighbors each point is connected to
    'Normalization' , 'bimarkov', ...
    'kEigenVecs'    , 500);                   % How to normalize diffusion matrix
G = GraphDiffusion(X, 0.05,graphDiffusion_opts);


figure;hold on;
scatter3(G.EigenVecs(:,2),G.EigenVecs(:,3),G.EigenVecs(:,4),30,G.EigenVecs(:,5),'filled');title('Diffusion map (2,3,4), color=next eigenvector');
scatter3(G.EigenVecs(:,3),G.EigenVecs(:,4),G.EigenVecs(:,5),30,G.EigenVecs(:,6),'filled');title('Diffusion map (3,4,5), color=next eigenvector');


figure;
gershdisc(G.W);
hold on;
d1 = eig(full(G.T));%,size(X,2),'SM');
Z=X;Z(2,:)=0;
Z(1,:)=d1;
plot(Z(1,:),Z(2,:),'y');



figure ;plot(X(1,:),X(2,:));


XY_fr=fruchterman_reingold_force_directed_layout(G.W);


figure;XY_KK=kamada_kawai_spring_layout(G.W,'tol',1e-6,'maxiter',500,'spring_constant',.1,...
    'progressive',0,'edge_length',1,'edge_weight','matrix');
[he,hv]=wgPlot(G.W,XY_fr,'vertexWeight',X(1,:),'vertexMetadata',X(2,:),'edgeColorMap',pink);



%%%%%%%%%Principle Curvature of Triangulated Mesh
load('testdata.mat');
[Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2]=patchcurvature(FV,true);
figure, title('Principal A');
p1=FV.vertices-2*Dir1; p2=FV.vertices+2*Dir1;       
plot3([p1(:,1) p2(:,1)]',[p1(:,2) p2(:,2)]',[p1(:,3) p2(:,3)]','g-');
axis equal; view(3) 
figure, title('Principal B');
p1=FV.vertices-2*Dir2; p2=FV.vertices+2*Dir2;       
plot3([p1(:,1) p2(:,1)]',[p1(:,2) p2(:,2)]',[p1(:,3) p2(:,3)]','r-');
axis equal; view(3)

load('testdata2.mat');
[Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2]=patchcurvature(FV);
figure,
subplot(2,2,1), title('Mean Curvature');
C=Cmean;
patch(FV,'FaceColor','interp','FaceVertexCData',C,'edgecolor','none');
axis equal; view(3)
subplot(2,2,2), title('Gaussian Curvature');
C=Cgaussian;
patch(FV,'FaceColor','interp','FaceVertexCData',C,'edgecolor','none');
axis equal; view(3)
subplot(2,2,3), title('Principal A');
p1=FV.vertices-2*Dir1; p2=FV.vertices+2*Dir1;       
plot3([p1(:,1) p2(:,1)]',[p1(:,2) p2(:,2)]',[p1(:,3) p2(:,3)]','g-');
axis equal; view(3) 
subplot(2,2,4), title('Principal B');
p1=FV.vertices-2*Dir2; p2=FV.vertices+2*Dir2;       
plot3([p1(:,1) p2(:,1)]',[p1(:,2) p2(:,2)]',[p1(:,3) p2(:,3)]','r-');
axis equal; view(3)





