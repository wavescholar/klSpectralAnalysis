function [ output_args ] = GrowClusteredGraph( Nodes, Clusters, InterclusterProbs,IntraClusterProbs)
Nodes=500;
N=Nodes;

x = randperm(N);
gs = N/Clusters;
Clusters=5;
InterclusterProbs=zeros(1,Nodes/Clusters);
mu=2;sigma=0.5;
ClusterProbs=normrnd(mu,sigma,Clusters,Nodes/Clusters);
InterClusterProbs=[.2,.2,.2,.2,.2]
%Make a struct with the cluster graphs G_C, and then mix to get  
stuct
for C=1:Clusters
probs=ClusterProbs(i,:);
G(C).
Cluster 
ClusterProbs(gs,gs) < InterClusterProbs(i);

end

return A;
end
