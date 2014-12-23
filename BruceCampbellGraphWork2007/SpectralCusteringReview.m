%  This is nice - revive this work


clear all;
N=3*  123;
x = (1:N);%randperm(N);
gs = N/3;
G1 = x(1:gs);
G2 = x(gs+1:2*gs);
G3= x(2*gs+1:3*gs);
p_G1 = 0.179;
p_G2 = 0.449;
p_G3 = 0.699;
p_Inbetween = 0.000;
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
A = triu(A,1);
A = A + A';
L2=laplacian(A);

[V D] = eigs(L2, 3, 'SA');
spy(A);
sigma=eig(L2);
hist(sigma,200);
plot(sort(V(:,2)), '.-');
plot(sort(V(:,3)), '.-');
Ic = adjacency2incidence(A);
Af = abs(Ic*Ic');
XY_fr=fruchterman_reingold_force_directed_layout(sparse(Af*1.0),'initial_temp',150);
figure ;wgPlot(Af,XY_fr,'vertexWeight',V(:,2),'edgecolormap',pink);


clear all;
Cox = [
	0.00  0.00  1.00
	0.00  0.50  0.00
	1.00  0.00  0.00
	0.00  0.75  0.75
	0.75  0.00  0.75
	0.75  0.75  0.00
	0.25  0.25  0.25
	0.75  0.25  0.25
	0.95  0.95  0.00
	0.25  0.25  0.75
	0.75  0.75  0.75
	0.00  1.00  0.00
	0.76  0.57  0.17
	0.54  0.63  0.22
	0.34  0.57  0.92
	1.00  0.10  0.60
	0.88  0.75  0.73
	0.10  0.49  0.47
	0.66  0.34  0.65
	0.99  0.41  0.23
];
for p=.05:.05:.5
    
    p_G1 = 0.25 + p;
    p_G2 = 0.25 + p;
    
    gaps = zeros(200,1);
    gcf = figure;
    %('Visible','off');
    hold on;
    clusterSizeStart = 128;
    clusterSizeDelta = 64;
    clusterSizeMax = 512;
    legendCount =1;
    for ClusterSize=clusterSizeStart:clusterSizeDelta:clusterSizeMax
        tic
        legendStrings{legendCount} = ['N = ',sprintf('%d',ClusterSize)];legendCount = legendCount+1;
        
        
        
        for i=1:200
            N=2*  ClusterSize;
            x = (1:N);%randperm(N);
            gs = N/2;
            G1 = x(1:gs);
            G2 = x(gs+1:2*gs);
            p_Inbetween = 0.000 + i/200;
            A(G1, G1) = rand(gs,gs) < p_G1;
            A(G2, G2) = rand(gs,gs) < p_G2;
            [n1,m1]=size( A(G1, G2));
            B_1=rand(gs, gs);
            A(G1, G2) = B_1 < p_Inbetween;
            A = triu(A,1);
            A = A + A';%spy(A);
            L2=laplacian(A);
            sigma=eig(L2);
            gap = sigma(2)-sigma(1);
            gaps(i)= gap;
            
            [V D] = eigs(L2, 3, 'SA');
            spy(A);
            sigma=eig(L2);
            hist(sigma,200);
            plot(sort(V(:,2)), '.-');
            plot(sort(V(:,3)), '.-');
            Ic = adjacency2incidence(A);
            Af = abs(Ic*Ic');
            XY_fr=fruchterman_reingold_force_directed_layout(sparse(Af*1.0),'initial_temp',150);
            figure ;wgPlot(Af,XY_fr,'vertexWeight',V(:,2),'edgecolormap',pink);
            
        end
        xVals = (1/200:1/200:200/200);
        %plot(xVals,gaps,'Color',[0+ .5*ClusterSize / clusterSizeMax ,0+ .5*ClusterSize / clusterSizeMax,0 + .5*ClusterSize / clusterSizeMax]);
        plot(xVals,gaps,'Color',Cox(legendCount,:));
        toc
    end
    legend(legendStrings,'Location','NW');
    title({['Spectral Gap For 2 Class Ensemble with varying inter class connectivity '];[' p(Class1)= ',sprintf('%f',p_G1),' p(Class2)= ',sprintf('%f',p_G2)]});
    xlabel('p(Class2  ~ Class1)');
    ylabel('\lambda_1 - \lambda_2');
    print(gcf,'-dpdf','-r600', ['K:\KL_A\klMatlabGraphTheory\SpectralCusteringReview_Figures\','SpectralGap_P1_',sprintf('%f',p_G1),'P2_',sprintf('%f',p_G2),'.pdf'] );
    close gcf;
end





























clear all;
Cox = [
	0.00  0.00  1.00
	0.00  0.50  0.00
	1.00  0.00  0.00
	0.00  0.75  0.75
	0.75  0.00  0.75
	0.75  0.75  0.00
	0.25  0.25  0.25
	0.75  0.25  0.25
	0.95  0.95  0.00
	0.25  0.25  0.75
	0.75  0.75  0.75
	0.00  1.00  0.00
	0.76  0.57  0.17
	0.54  0.63  0.22
	0.34  0.57  0.92
	1.00  0.10  0.60
	0.88  0.75  0.73
	0.10  0.49  0.47
	0.66  0.34  0.65
	0.99  0.41  0.23
];
for p=.05:.05:.5
    
    p_G1 = 0.25 + p;
    p_G2 = 0.25 + p;
    
    gaps12 = zeros(200,1);
    gaps23 = zeros(200,1);
    gcf = figure;
    %('Visible','off');
    hold on;
    clusterSizeStart = 128;
    clusterSizeDelta = 64;
    clusterSizeMax = 512;
    legendCount =1;
    for ClusterSize=clusterSizeStart:clusterSizeDelta:clusterSizeMax
        tic
        legendStrings{legendCount} = ['N = ',sprintf('%d',ClusterSize)];legendCount = legendCount+1;

       for i=1:200
            N=3*  ClusterSize;
            x = (1:N);%randperm(N);
            gs = N/3;
            G1 = x(1:gs);
            G2 = x(gs+1:2*gs);
            G3= x(2*gs+1:3*gs);
            p_G1 = 0.079;
            p_G2 = 0.0449;
            p_G3 = 0.0699;
            p_Inbetween = 0.000+ i/10000;
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
            A = triu(A,1);
            A = A + A';

            L2=laplacian(A);
            sigma=eig(L2);
            gap = sigma(2)-sigma(1);
            gaps12(i)= gap;
            
            gap = sigma(3) - sigma(2);
            gaps23(i) = gap;
            
            [V D] = eigs(L2, 3, 'SA');
            figure;spy(A);
%             sigma=eig(L2);
%             hist(sigma,200);
           plot(sort(V(:,2)), '.-');
           plot(sort(V(:,3)), '.-');
            Ic = adjacency2incidence(A);
            Af = abs(Ic*Ic');
            XY_fr=fruchterman_reingold_force_directed_layout(sparse(Af*1.0),'initial_temp',150);
            figure ;wgPlot(Af,XY_fr,'vertexWeight',V(:,2),'edgecolormap',pink);
            
        end
        xVals = (1/200:1/200:200/200);
        %plot(xVals,gaps,'Color',[0+ .5*ClusterSize / clusterSizeMax ,0+ .5*ClusterSize / clusterSizeMax,0 + .5*ClusterSize / clusterSizeMax]);
        set(0, 'DefaultAxesLineStyleOrder','-|-.')

        plot(xVals,gaps12,'Color',Cox(legendCount,:));hold on;
        plot(xVals,gaps23,'Color',Cox(legendCount+1,:));
        toc
    end
    legend(legendStrings,'Location','NW');
    title({['Spectral Gap For 2 Class Ensemble with varying inter class connectivity '];[' p(Class1)= ',sprintf('%f',p_G1),' p(Class2)= ',sprintf('%f',p_G2)]});
    xlabel('p(Class2  ~ Class1)');
    ylabel('\lambda_1 - \lambda_2');
    print(gcf,'-dpdf','-r600', ['K:\KL_A\klMatlabGraphTheory\SpectralCusteringReview_Figures\','SpectralGap_P1_',sprintf('%f',p_G1),'P2_',sprintf('%f',p_G2),'.pdf'] );
    close gcf;
end