function [costs] = mdijkstra(A,C)
%
%A=square matrix (either adjecancy or cost)
%
%if C=1 then A=adjecancy matrix
%    where, element(i,j)=1 when vertex v is directly connected with j
%    else (i,j)=0
%
%if C=2 then A=cost matrix 
%    where, element (i,j) represents positive integer representing cost
%    between vertex i and j
%
% Output: [costs]: calculated cost matrix
% Developed by: Bharat Patel
% Release date: 03/28/2009

n = length(A);
costs=single(A);
costs(costs==0)=Inf;
for k = 1:n
    disp(sprintf('%d/13-%d/%d',C,k,n));
    w_col=single(costs(:,k));
    if C==1
        mx=max(w_col(find(w_col~=Inf)));
        pmx=0;
        while(mx~=pmx)
            cols=find(w_col==mx);
            tmp=min(costs(:,cols),[],2);
            tmp=tmp+mx;
            tmp1=[w_col tmp];
            w_col=min(tmp1,[],2);
            pmx=mx;
            mx=max(w_col(find(w_col~=Inf)));
        end
        costs(:,k)=w_col;
    elseif C==2
        m1=(w_col(find(w_col)));
        m=sort(unique(m1),'ascend');
        for j=1:length(m)
            mx=m(j);
            cols=find(w_col==mx);
            tmp=min(costs(:,cols),[],2);
            tmp=tmp+mx;
            tmp1=[w_col tmp];
            w_col=min(tmp1,[],2);
        end
        costs(:,k)=w_col;
        costs(k,:)=w_col';
    end
end
for k=1:n
    costs(k,k)=0;
end



