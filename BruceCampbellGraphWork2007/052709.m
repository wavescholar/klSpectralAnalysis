% Three cluseters with a cut between each
clear all;
N=3*    129;
% In this segment, we'll plant a partition in a graph, and then use the second smallest eigenvector to find it.
% As always, the first step is to generate our dataset. In this example, we'll be a little more ambitious and use a larger number of vertices. 
x = (1:N);%randperm(N);
gs = N/3;
G1 = x(1:gs);
G2 = x(gs+1:2*gs);
G3= x(2*gs+1:3*gs);
p_G1 = 0.13;
p_G2 = 0.23;
p_G3 = 0.09;
p_Inbetween = 0.008;
A(G1, G1) = rand(gs,gs) < p_G1;
A(G2, G2) = rand(gs,gs) < p_G2;
A(G3, G3) = rand(gs,gs) < p_G3;

[n1,m1]=size( A(G1, G2));
[n2,m2]=size( A(G2, G3));
[n3,m3]=size( A(G3, G1));

B_1=rand(gs, gs);B_2=rand(gs, gs);B_3=rand(gs, gs);
 A(G1, G2) = B_1 < p_Inbetween;
 A(G2, G3) = B_2 < p_Inbetween/3;
 A(G1, G3) = B_3 < p_Inbetween/5;

A = triu(A,1);
A = A + A';
spy(A);

%L = del2(A);
L2=laplacian(A);%

Ic = adjacency2incidence(A);
%L3 = Ic*Ic';%compute_combinatorial_laplacian(A);

% spy(L);
% spy(L2);
% spy(L3-L2);

[V D] = eigs(L2, 3, 'SA');sigma=eig(A);hist(sigma,50);
D(2,2)
plot(V(:,2), '.-');
plot(sort(V(:,2)), '.-');
plot(V(:,3), '.-');
plot(sort(V(:,3)), '.-');

spy(A);%draw_dot(A);
[ignore p] = sort(V(:,2));
spy(A(p,p));

XY_ga = gursoy_atun_layout(sparse(A*1.0),'topology','circle');
XY_kk = kamada_kawai_spring_layout(sparse(A*1.0));
XY_fr=fruchterman_reingold_force_directed_layout(sparse(A*1.0));
% wgplot(GEV_l,X);
% AB=sparse(double(Net));
% XB = gursoy_atun_layout(AB,'topology','heart');gplot(AB,X,'.-');

wgplot(A,XY_ga);
wgplot(A,XY_kk);
wgplot(A,XY_fr); %This is the best format for cut display
gplot(A,XY_fr);
WDEFAULT(gcf)
[row_clust_idx, col_clust_idx,y_index,x_index]=SpectralCoClustering(A,3);%,display,names)


[ w_fdla ] = fdla(sparse(A));
[ w_fmmc ] = fmmc(Ic);
[ w_md   ] = max_deg(Ic);
[ w_bc   ] = best_const(A);
[ w_mh   ] = mh(A);

Afmmc = abs(Ic*diag(w_fmmc)*Ic');
Amd = abs(Ic*diag(w_md)*Ic');
An = (An - diag(diag(An))) > 0.001;
XY_fr=fruchterman_reingold_force_directed_layout(sparse(Afmmc*1.0));
wgplot(Afmmc,XY_fr,'vertexWeight',V(:,3));



XY_fr=fruchterman_reingold_force_directed_layout(sparse(Afmmc*1.0));
wgplot(Afmmc,XY_fr,'vertexWeight',V(:,3));


XY_fr=fruchterman_reingold_force_directed_layout(sparse(A));
wgplot(A,XY_fr,'vertexWeight',w_fdla);

cm=[w_fmmc,w_fmmc,w_fmmc];
tic;
[he,hv]=wgPlot(Afmmc,XY_fr,'vertexWeight',V(:,3),'edgecolormap',pink);
toc


%Generate sample from \beta Hermite ensemble
n=N;m=N;
beta=1;
d=sqrt(chi2rnd(beta*[n:-1:1]))';
H=spdiags(d,1,n,n)+spdiags(randn(n,1),0,n,n);
H=(H+H')/2;
[L Q]=eigs(H);
hist(L);v=diag(H);


gplot(An,XY_fr,'b-'); hold on;
h = findobj(gca,'Type','line');
set(h,'LineWidth',1.5)
gplot(A,XY_fr,'b:');
plot(xy(:,1), xy(:,2), 'ko','LineWidth',4, 'MarkerSize',4);
axis([0.05 1.1 -0.1 0.95]);
title('Subgraph with positive transition prob.')
hold off;


BB = triu(B_1,1);
BB = BB + BB';
 [W,S] = FMMC(BB) %gives a vector of the fastest mixing Markov chain
% edge weights for a graph described by the incidence matrix A (n x m).
% Here n is the number of nodes and m is the number of edges in the graph;
% each column of A has exactly one +1 and one -1.
%
% The FMMC edge weights are given the SDP:
%
%   minimize    s
%   subject to  -s*I <= I - L - (1/n)11' <= s*I
%               w >= 0,  diag(L) <= 1
%


%Maeshpart Pkg
specpart(A,XY_fr);
cutsize(A,ans);

geopart(A,XY_fr,2);
cutsize(A,ans);

metispart(sparse(A),sparse(XY_fr));
cutsize(A,ans)

gspart(A,XY_fr);
disp('cutsize(Tapir,ans)');
cutsize(A,ans)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MESHPART

%%%%%%%%%%%%%%%%%%%%%%%%QUANTUM CLUSTERING - get the NN toolbox to run this
% % function [mapping] = qcClustering(sigma,steps,data)
% % performs quantum clustering + gradient descent
% % Input -
% % sigma - as in the algorithm
% % steps - number of steps
% % data - dims*samples matrix
% % doRecurse - if 1, then call the methods again - if not use the
% % bestCluster
% % numOfClusters - number of clusters
% % Output -
% % mapping results of the algorithm
%  [mapping] = qcClustering(2,20,XY_fr,1,1,3)









C_coeff=clustering_coefficients((sparse(A)));full(A)
hist(C_coeff,100);C_mu=mean(C_coeff);

Net = SFNG(N+20, 5, A);
C_coeff=clustering_coefficients(sparse(double(Net)));
hist(C_coeff,100);C_mu=mean(C_coeff);

[x y]=draw_dot(Net);
gplot(Net, [x' y'], '.-');












clear all;
N=150;
% In this segment, we'll plant a partition in a graph, and then use the second smallest eigenvector to find it.
% As always, the first step is to generate our dataset. In this example, we'll be a little more ambitious and use a larger number of vertices. 
x = randperm(N);
gs = N/2;
G1 = x(1:gs);
G2 = x(gs+1:end);
% Decide on the probabilities of edges within each group and between the two groups. 
% Because we are planting a partition, the probabilities of edges between the groups should be much 
% lower than the probability of edges within each group. Suppose that group 1 is a little more tightly connected than group 2.
p_G1 = 0.23;
p_G2 = 0.23;
p_Inbetween = 0.04;
A(G1, G1) = rand(gs,gs) < p_G1;
A(G2, G2) = rand(N-gs,N-gs) < p_G2;
A(G1, G2) = rand(gs, N-gs) < p_Inbetween;
At = triu(A,1);

A = At -At';
Ai=incidence(A);

Aa=incidence2adjacency(A);

L = del2(A);
L2=laplacian(A);%

Ic = adjacency2incidence(A);
L3 = Ic*Ic';%compute_combinatorial_laplacian(A);




[V D] = eigs(L2, 2, 'SA');
D(2,2);%ans = 46.7158
plot(V(:,2), '.-');
plot(sort(V(:,2)), '.-');
[ignore p] = sort(V(:,2));
XY_fr=fruchterman_reingold_force_directed_layout(sparse(At+At'*1.0));

wgplot(A,XY_fr); %This is the best format for cut display

xy=XY_fr;
%%%%%%%%%%%%%%%%%%%%%Meshpart Pkg
specpart(A,XY_fr);
cutsize(A,ans);

geopart(A,XY_fr,2);
cutsize(A,ans);

metispart(sparse(A),sparse(XY_fr));
cutsize(A,ans)

gspart(A,XY_fr);
disp('cutsize(Tapir,ans)');
cutsize(A,ans)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MESHPART


[ w_fdla ] = fdla(sparse(A));
[ w_fmmc ] = fmmc(Aa);
[ w_md   ] = max_deg(A);
[ w_bc   ] = best_const(A);
[ w_mh   ] = mh(A);


[he,hv]=wgPlot(A,XY_fr,'vertexWeight',w_fdla);
[he,hv]=wgPlot(A,XY_fr,'vertexWeight',w_fmmc);
[he,hv]=wgPlot(A,XY_fr,'vertexWeight',w_md);
[he,hv]=wgPlot(A,XY_fr,'vertexWeight',w_bc);
[he,hv]=wgPlot(A,XY_fr,'vertexWeight',w_mh);

[he,hv]=wgPlot(A,XY_fr,'vertexWeight',w_fmmc,'vertexMetadata',w_fmmc*5,'vertexScale',200,'edgeColorMap',pink);
%bbcrevisist lint fix dunno
%rand(length(A),1),'edgeColorMap';pink));











clear all;
memory;%clear, pack, whos, inmem, save, load, mlock, munlock
N=5*132;
% In this segment, we'll plant a partition in a graph, and then use the second smallest eigenvector to find it.
% As always, the first step is to generate our dataset. In this example, we'll be a little more ambitious and use a larger number of vertices. 
x = (1:N);%randperm(N);
gs = N/5;

G1 = x(1:gs);
G2 = x(gs+1:2*gs);
G3 = x(2*gs+1:3*gs);
G4 = x(3*gs+1:4*gs);
G5 = x(4*gs+1:5*gs);
%V= exp(1)./ log(1:5*N);plot(V);
V= (1)./ log(1:5*N);plot(V);
p_G1 = 0.53*V(1*N);
p_G2 = 0.43*V(1*N);
p_G3 = 0.39*V(1*N);
p_G4 = 0.29*V(1*N);
p_G5 = 0.19*V(1*N);

p_Inbetween_12 = 0.016*V(2*N);
p_Inbetween_13 = 0.010*V(2*N);
p_Inbetween_14 = 0.0108*V(2*N);
p_Inbetween_15 = 0.0106*V(2*N);

p_Inbetween_23 = 0.018*V(2*N);
p_Inbetween_24 = 0.0106*V(2*N);
p_Inbetween_25 = 0.0104*V(2*N);

p_Inbetween_34 = 0.018*V(2*N);
p_Inbetween_35 = 0.0106*V(2*N);

p_Inbetween_45 = 0.016;

A(G1, G1) = rand(gs,gs) < p_G1;
A(G2, G2) = rand(gs,gs) < p_G2;
A(G3, G3) = rand(gs,gs) < p_G3;
A(G4, G4) = rand(gs,gs) < p_G4;
A(G5, G5) = rand(gs,gs) < p_G5;

B_1=rand(gs, gs);
B_2=rand(gs, gs);
B_3=rand(gs, gs);
B_4=rand(gs, gs);
B_5=rand(gs, gs);

A(G1, G2) = B_1 < p_Inbetween_12;
A(G1, G3) = B_1 < p_Inbetween_13;
A(G1, G4) = B_1 < p_Inbetween_14;
A(G1, G5) = B_1 < p_Inbetween_15;


A(G2, G3) = B_2 < p_Inbetween_23;
A(G2, G4) = B_2 < p_Inbetween_24;
A(G2, G5) = B_2 < p_Inbetween_25;


A(G3, G4) = B_3 < p_Inbetween_34;
A(G3, G5) = B_3 < p_Inbetween_35;

A(G4, G5) = B_4 < p_Inbetween_45;

A = triu(A,1);
A = A + A';
spy(A);

%L = del2(A);
%L2=laplacian(A);%

Ic = adjacency2incidence(A);
L3 = Ic*Ic';%compute_combinatorial_laplacian(A);
[V D] = eigs(L3, 5, 'SA');
sigma=eig(A);
hist(sigma,50);
D(2,2)

plot(V(:,2), '.-');
plot(sort(V(:,2)), '.-');

plot(V(:,3), '.-');
plot(sort(V(:,3)), '.-');

plot(V(:,4), '.-');
plot(sort(V(:,4)), '.-');

plot(V(:,5), '.-');
plot(sort(V(:,5)), '.-');


spy(A);%draw_dot(A);
[ignore p] = sort(V(:,2));
spy(A(p,p));

XY_fr=fruchterman_reingold_force_directed_layout(sparse(A*1.0));


wgplot(A,XY_fr); %This is the best format for cut display
gplot(A,XY_fr);

memory

C_coeff=clustering_coefficients(sparse(A));
plot(hist(C_coeff,100));

[ w_fdla ] = fdla(sparse(A));
[ w_fmmc ] = fmmc(Ic);
[ w_md   ] = max_deg(Ic);
[ w_bc   ] = best_const(A);
[ w_mh   ] = mh(A);

Afmmc = abs(Ic*diag(w_fmmc)*Ic');
%Amd = abs(Ic*diag(w_md)*Ic');
% Comment this out to keep weighted adjacency matrix

An = (An - diag(diag(An))) > 0.001;
XY_fr=fruchterman_reingold_force_directed_layout(sparse(Afmmc*1.0));
wgplot(Afmmc,XY_fr,'vertexWeight',V(:,3));


%Waxman random network topology generator
%The random network topology generator introduced in Waxman (1988) is a geographic model for the growth of a computer network. In this model the nodes of the network are uniformly distributed in the plane and edges are added according to probabilities that depend on the distances between the nodes. The probability to have an edge between nodes u and v is given by 
%P(u, v) = ae-d/(bL)  
lambda=2,alpha=1, beta=.5, domain=[-1 1 -1 1] .* 3;
[adj_matr, nd_coord] = waxtop(lambda,alpha, beta, domain); 
C_coeff=clustering_coefficients(sparse(A));
plot(hist(C_coeff,100));

% Inputs: 
%   lambda - intensity of the Poisson process
%   alpha - maximal link probability
%   beta - parameter to control length of the edges. Increased <beta>
%     yields a larger ratio of long edges to short edges
%   domain - bounds for the region. A 4-dimensional vector in
%     the form [x_min x_max y_min y_max]


% %Galton-Watson 
ngen=7;
p=[0 1 1 1 2 3 4 3 4 3 4 6 7 8 10 10 11 11 11 12 15 16 18 17 17 17 21 23 24 24 22 25 22 25 22 25 27 28 34 31 31 29 30 36 29 30 36 29 30 36];
[parents] = branch(ngen, p);

% There was concern amongst the Victorians that aristocratic surnames were becoming extinct. Galton originally posed the question regarding the probability of such an event in the Educational Times of 1873, and the Reverend Henry William Watson replied with a solution. Together, they then wrote an 1874 paper entitled On the probability of extinction of families. Galton and Watson appear to have derived their process independently of the earlier work by I. J. Bienaymé; see Heyde and Seneta 1977. For a detailed history see Kendall (1966 and 1975).
% 
% 
% [edit] Concepts
% Assume, as was taken for granted in Galton's time, that surnames are passed on to all male children by their father. Suppose the number of a man's sons to be a random variable distributed on the set { 0, 1, 2, 3, ...}. Further suppose the numbers of different men's sons to be independent random variables, all having the same distribution.
% 
% Then the simplest substantial mathematical conclusion is that if the average number of a man's sons is 1 or less, then their surname will surely die out, and if it is more than 1, then there is more than zero probability that it will survive forever.
% 
% Modern applications include the survival probabilities for a new mutant gene, or the initiation of a nuclear chain reaction, or the dynamics of disease outbreaks in their first generations of spread, or the chances of extinction of small population of organisms; as well as explaining (perhaps closest to Galton's original interest) why only a handful of males in the deep past of humanity now have any surviving male-line descendants, reflected in a rather small number of distinctive human Y-chromosome DNA haplogroups.
% 
% A corollary of high extinction probabilities is that if a lineage has survived, it is likely to have experienced, purely by chance, an unusually high growth rate in its early generations at least when compared to the rest of the population.
% 
% 
% [edit] Mathematical definition
% A Galton-Watson process is a stochastic process {Xn} which evolves according to the recurrence formula X0 = 1 and
% 
%  
% where for each n,  is a sequence of IID natural number-valued random variables. The extinction probability (i.e. the probability of final extinction)is given by
% 
%  
% This is clearly equal to zero if each member of the population has exactly one descendent. Excluding this case (usually called the trivial case) there exists a simple necessary and sufficient condition, which is:
% 
% 
% [edit] Extinction criterion for Galton-Watson process
% In the non-trivial case the probability of final extinction is equal to one if E{?1} ? 1 and strictly less than one if E{?1} > 1.
% 
% The process can be treated analytically using the method of probability generating functions.
% 
% If the number of children ? j at each node follows a Poisson distribution, a particularly simple recurrence can be found for the total extinction probability xn for a process starting with a single individual at time n = 0:
% 
%  
% giving the curves plotted above.
% 
% 
% [edit] Bisexual Galton–Watson process
% In the (classical) Galton–Watson process defined above, only men count, that is, the reproduction can be understood as being asexual. The more natural corresponding version for (bi)sexual reproduction is the so-called 'bisexual Galton–Watson process', where only couples can reproduce. In this process, each child is supposed to be male or female, independently of each other, with a specified probability, and a so-called 'mating function' determines how many couples will form in a given generation. As before, reproduction of different couples are considered to be independent of each other. Now the analogue of the trivial case corresponds to the case of each male and female reproducing in exactly one couple, having one male and one female descendent, and that the mating function takes the value of the minimum of the number of males and females (which are then the same from the next generation onwards).
% 
% Since the total reproduction within a generation depends now strongly on the mating function, there exists in general no simple necessary and sufficient condition for final extinction as it is the case in the classical Galton–Watson process. However, excluding the non-trivial case, the concept of the averaged reproduction mean (Bruss (1984)) allows for a general sufficient condition for final extinction, treated in the next section.
% 
% 
% [edit] Extinction criterion (bisexual Galton-Watson process)
% If in the non-trivial case the averaged reproduction mean per couple stays bounded over all generations and will not exceed 1 for a sufficiently large population size, then the probability of final extinction is always 1.
% 
% 
% [edit] Examples
% Countries that have used family names for many generations exhibit the Galton–Watson process in their low number of surviving family names:
% 
% Korean names are the most striking example, with 250 family names, and 45% of the population sharing three family names 
% Chinese names are similar, with 22% of the population sharing three family names (numbering close to 300 million people), and the top 200 names covering 96% of the population. 
% By contrast:
% 
% Dutch names have only included a family name since the Napoleonic Wars in the early 19th century, and there are over 68,000 Dutch family names. 
% Thai names have only included a family name since 1920, and only a single family can use a given family name, hence there are a great number of Thai names. Further, Thai people change their family names with some frequency, complicating the analysis. 
% 
% [edit] See also
% Branching process 
% 
% [edit] References
% F. Thomas Bruss (1984). "A Note on Extinction Criteria for Bisexual Galton–Watson Processes". Journal of Applied Probability 21: 915–919. 
% C C Heyde and E Seneta (1977). I.J. Bienayme: Statistical Theory Anticipated. Berlin, Germany. 
% D G Kendall (1966). Journal of the London Mathematical Society 41: 385–406 
% D G Kendall (1975). Bulletin of the London Mathematical Society 7: 225–253 
% H W Watson and Francis Galton, "On the Probability of the Extinction of Families", Journal of the Anthropological Institute of Great Britain, volume 4, pages 138–144, 1875. 
% 
% [edit] External links
% The original Galton–Watson paper: On the Probability of the Extinction of Families 
% "Survival of a Single Mutant" by Peter M. Lee of the University of York 
% Retrieved from "http://en.wikipedia.org/wiki/Galton%E2%80%93Watson_process"
