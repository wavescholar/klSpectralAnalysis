function [Yrows nclasses] = gingcca_symmetric(A,tol)
% segments the affinity matrix A using cca
% define the neighborhood of a(j,k) as da(j,k) = {row j and column k} =
% {a(j,:) and a(:,k)}
% then the set of neighbors that are connected to a(j,k) is: 
% c(a(j,k)) = {a(u,v) in da(j,k) : a(u,v) > tol} 

% if want actual labeled matrix (Y), use gingcca_speedup.m (now slower than
% this)

% assumes A is positive!!! only classifies points (j,k) with A(j,k) > tol

% algorithm from
% https://engineering.purdue.edu/~bouman/ece637/notes/pdb/ConnectComp.pdf 
% and c(s) designed for affinity matrix by me
disp('assumes A is symmetric!!!');

if nargin < 2
    tol = 0;  % connected if A(i,j) > 0.
end

% start with a vector of non-classified points
[nj nk] = size(A);
Ycols = zeros(1,nk,'int32'); % preliminary labels for each column
Yrows = zeros(nj,1,'int32'); % preliminary labels for each row
%Y = sparse(nj,nk); % vector of class labels (Y = 0 --> not yet classified)
nclasses = 0; % no classes have been formed yet

ncl = 0; % keep track of the number of columns labeled, to pass the time

% start by giving all non-zero (> tol) entries in each column a label
Cols = cell(nk,1);
%Rows = cell(nj,1);
for k = 1:nk
    if mod(k,10000) == 1
        disp([num2str(k),' : ',num2str(nk)]);
    end
    Cols{k} = int32(find(A(:,k) > tol)); % all the non-zero row indices in column k
    %Rows{k} = Cols{k}';
end
disp('found all non-zeros in columns - also in rows, bc assuming symmetric!');

% Rows = cell(nj,1);
% for j = 1:nj
%     if mod(j,10000) == 1
%         disp([num2str(j),' : ',num2str(nj)]);
%     end
%     % take advantage of symmetry of A!!
%     Rows{j} = Cols{j}'; % find(A(j,:) > tol); 
% end
% disp('found all non-zeros in rows');

maxbend = 0;

disp('starting to connect rows and columns');
% now start connecting columns and rows
for k = 1:nk
    if Ycols(k) == 0 % if this column hasn't been labeled yet
        B = zeros(nj,3,'int32'); % = cell(nj*nk,1);
        lengthB = nj;
        b1 = 1;
        bend = 0;
        rc = 2; %column
        for c = 1:length(Cols{k})
            bend = bend + 1; % add an entry to B
            B(bend,:) = [Cols{k}(c) k rc]; % B{bend} = [Cols{k}(c) k rc]
        end
        nB = bend - b1 + 1;
        
        if nB > 0
            nclasses = nclasses + 1;
        end
        
        while nB > 0
            
            if bend > lengthB
                disp('storage matrix B too small: doubling in length');
                newB = zeros(2*size(B,1),3,'int32');
                newB(1:size(B,1),:) = B;
                B = newB;
                lengthB = size(B,1);
                clear newB;
            end
            
            s = B(b1,:); % = B{b1};
            u = s(1);
            v = s(2);
            rc = s(3);
            B(b1,:) = [0 0 0];
            %Y(u,v) = nclasses; % label this point
            b1 = b1 + 1; % leave out s from range of current items in B
            % add all entries in s's row or column to B
            if rc == 2 % column, so get row data
                if Yrows(u) == 0 % if this row hasn't been  labeled yet
                    for c = 1:length(Cols{u})
                        bend = bend + 1;
                        B(bend,:) = [u Cols{u}(c) 1]; % B{bend} = [u Cols{u}(c) 1];
                    end
                    Yrows(u) = nclasses;
                end
            end
            if rc == 1 % row, so get column data
                if Ycols(v) == 0
                    for c = 1:length(Cols{v})
                        bend = bend + 1;
                        B(bend,:) = [Cols{v}(c) v 2]; %  B{bend} = [Cols{v}(c) v 2];
                    end
                    Ycols(v) = nclasses;
                    ncl = ncl + 1;
                    if mod(ncl,100) == 1
                        disp([num2str(ncl),' : ',num2str(nk),' columns labeled']);
                    end
                    
                end
            end
            nB = bend - b1 + 1; % get size of current indices in B
            
        end % end while nB > 0
        if bend > maxbend
        maxbend = bend;
        end
    end
end

maxbend % just to figure out how large to preallocate the matrix B