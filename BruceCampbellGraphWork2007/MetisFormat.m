function y=MetisFormat(A,filename)

% This fn convers an Adjacency matrix to metis format

% The first line contains either two (n, m), three (n, m, fmt), or four (n, m, fmt, ncon) integers. The first two integers
% (n, m) are the number of vertices and the number of edges, respectively. Note that in determining the number of edges
% m, an edge between any pair of vertices v and u is counted only once and not twice (i.e., we do not count the edge
% .v; u/ separately from .u; v/). For example, the graph in Figure 8 contains 11 vertices. The third integer (fmt) is used
% to specify whether or not the graph has weights associated with its vertices, its edges, or both. Table 2 describes the
% possible values of fmt and their meaning. Note that if the graph is unweighted (i.e., all vertices and edges have the
% same weight), then the fmt parameter can be omitted. Finally, the fourth integer (ncon) is used to specify the number
% of weights associated with each vertex of the graph. The value of this parameter determines whether or not METIS will
% use the multi-constraint partitioning algorithms described in Section 3. If the vertices of the graph have no weights or
% only a single weight, then the ncon parameter can be omitted. However, if ncon is greater than 0, then the file should
% contain the required vertex weights and the fmt parameter should be set appropriately (i.e., it should be set to either 10
% or 11).
% fmt Meaning
% 0 The graph has no weights associated with either the edges or the vertices
% 1 The graph has weights associated with the edges
% 10 The graph has weights associated with the vertices
% 11 The graph has weights associated with both the edges & vertices

%input the Adjacency matrix


%Test - see metis doc
% 7 11
% 5 3 2
% 1 3 4
% 5 4 2 1
% 2 3 6 7
% 1 3 6
% 5 4 7
% 6 4

    

filename='C:\KL\output\UnitTestOutput\graphTheory\test.graph'

% A=zeros(7,7); A test matrix from the metis documentation
% A(1, 5)=1; A(1,3)=1; A(1,2)=1;
% A(2,1)=1; A(2,3)=1; A(2,4)=1;
% A(3,5)=1; A(3,4)=1; A(3,2)=1 ;A(3,1)=1;
% A(4,2)=1; A(4,3)=1;A(4,6)=1; A(4,7)=1;
% A(5,1)=1 ;A(5,3)=1; A(5,6)=1;
% A(6,5)=1; A(6,4)=1; A(6,7)=1;
% A(7,6)=1 ;A(7,4)=1;

%All is this is to get the edge count
xy= kamada_kawai_spring_layout(sparse(A*1.0));
[V,E] = axy2ve(A,xy);
[l,p]=size(E);
fid=fopen(filename,'wt+');


k=0;[n,m]=size(A);
for i=1:n
   
    for j=1:m
        if A(i,j)==1
            k=k+1;
        end
    end
end
    
%Write Metis Header Info
fprintf(fid, '%u %u \n', n,l*2);
%Loop over each vertex i and make a list of the vertices connecting to i
for i=1:n
    k=0;
    for j=1:m
        if A(i,j)==1
            k=k+1;
        end
    end
    v=zeros(k,1);
    k=1;
    for j=1:m
        if A(i,j)==1
            v(k)=j;
            fprintf(fid, '%d  ', v(k));
            k=k+1;
            
            
        end
    end
    fprintf(fid, '\n');
    %count = fwrite(fid, v','double');
end
fclose(fid);