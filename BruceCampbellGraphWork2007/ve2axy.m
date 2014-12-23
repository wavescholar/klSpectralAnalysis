function [A,xy] = ve2axy(V,E)
% VE2AXY Convert Graph of Vertices and Edges to Adjacency Matrix and XY Points
%
% Inputs:
%     V   is a Nx2 (or Nx3) matrix of x,y,(z) coordinates
%     E   is a Px2 matrix containing a list of edge connections
%
% Outputs:
%     A   is a NxN adjacency matrix, where A(I,J) is nonzero
%           if and only if an edge connects point I to point J
%     xy  is a Nx2 (or Nx3) matrix of x,y,(z) coordinates (equivalent to V)
%
% Example:
%     n = 10;
%     V = 10*rand(n,2)
%     E = ceil(n*rand(3*n,2))
%     [A,xy] = ve2axy(V,E)
%     spy(A);
%
% Example:
%     n = 2e4;
%     V = 10*rand(n,2);
%     E = ceil(n*rand(n,2));
%     [A,xy] = ve2axy(V,E);
%     spy(A);
%
% Web Resources:
%   <a href="http://en.wikipedia.org/wiki/Graph_%28mathematics%29">Graphs</a>
%   <a href="http://en.wikipedia.org/wiki/Adjacency_matrix">Adjacency Matrix</a>
%
% See also: axy2ve, gplot, gplotd, gplotdc, distmat
%
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
% Release: 1.0
% Release Date: 5/22/08

error(nargchk(2,2,nargin));

n = size(V,1);
E = E(:,1:2);
if (min(E(:)) < 1) || (max(E(:)) > n)
    error('Invalid Edge List: Unexpected value not in the range 1 to %d.',n);
end

xy = V;
I = E(:,1);
J = E(:,2);
IJ = I + n*(J-1);

try
    A = zeros(n);
    A(IJ) = 1;
catch me
    try
        disp('Warning: Memory Limit Exceeded - Attempting alternate method.');
        A = spaxy();
    catch me2
        disp(me);
        rethrow(me2);
    end
end

    function spA = spaxy()
        % Create Sparse Adjacency Matrix and Process Edges in Blocks
        spA = sparse(n,n);
        m = length(IJ);
        block_size = max(1e3,min(1e4,1e2*floor(m/1e4)));
        N = floor(m/block_size);
        disp('Processing Edges in Blocks ...');
        for k = 1:N
            spA(IJ(block_size*(k-1)+1:block_size*k)) = 1;
        end
        spA(IJ(block_size*N+1:m)) = 1;
        disp('Warning: Adjacency matrix is SPARSE');
    end
end
