

clear all;
N=3*256;
% In this segment, we'll plant a partition in a graph, and then use the second smallest eigenvector to find it.
% As always, the first step is to generate our dataset. In this example, we'll be a little more ambitious and use a larger number of vertices. 
x = randperm(N);%;(1:N)
gs = N/3;
G1 = x(1:gs);
G2 = x(gs+1:2*gs);
G3= x(2*gs+1:3*gs);
p_G1 = 0.143;
p_G2 = 0.23;
p_G3 = 0.43;
p_Inbetween = 0.03;
A(G1, G1) = rand(gs,gs) < p_G1;
A(G2, G2) = rand(gs,gs) < p_G2;
A(G3, G3) = rand(gs,gs) < p_G3;

A=double(A);
[n1,m1]=size( A(G1, G2));
[n2,m2]=size( A(G2, G3));
[n3,m3]=size( A(G3, G1));

B_1=rand(gs, gs);B_2=rand(gs, gs);B_3=rand(gs, gs);
B_1=B_1 < p_Inbetween;B_2=B_2 < p_Inbetween/3;B_3=B_3 < p_Inbetween/5;
 A(G1, G2) = B_1;% B_1 < p_Inbetween;
 A(G2, G3) = B_2;%B_2 < p_Inbetween/3;
 A(G1, G3) = B_3;%B_3 < p_Inbetween/5;

A = triu(A,1);
A = A + A';spy(A);
L = del2(A);spy(L);
[V D] = eigs(L, 2, 'SA');sigma=eig(A);hist(sigma,50);
D(2,2)
plot(V(:,2), '.-');
plot(sort(V(:,2)), '.-');spy(A);%draw_dot(A);
[ignore p] = sort(V(:,2));
spy(A(p,p));%draw_dot(A(p,p));
C_coeff=clustering_coefficients(sparse(A));full(A)
hist(C_coeff,100);C_mu=mean(C_coeff);

Net = SFNG(N+20, 5, A);
C_coeff=clustering_coefficients(sparse(double(Net)));
hist(C_coeff,100);C_mu=mean(C_coeff);draw_dot(Net);


    sA=(sparse((A)));
    X1 = gursoy_atun_layout((sA),'topology','heart');gplot(sA,X1,'k');
    





