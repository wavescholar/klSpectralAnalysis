function [V,E] = axy2ve(A,xy)
% AXY2VE Convert Graph of Adjacency Matrix and XY Points to Vertices and Edges
%
% Inputs:
%     A   is a NxN adjacency matrix, where A(I,J) is nonzero
%           if and only if an edge connects point I to point J
%     xy  is a Nx2 (or Nx3) matrix of x,y,(z) coordinates (equivalent to V)
%
% Outputs:
%     V   is a Nx2 (or Nx3) matrix of x,y,(z) coordinates
%     E   is a Px2 matrix containing a list of edge connections
%
% Example:
%     n = 10;
%     xy = 10*rand(n,2)
%     A = round(rand(n))
%     spy(A);
%     [V,E] = axy2ve(A,xy)
%
% Example:
%     n = 2e4;
%     xy = 10*rand(n,2);
%     A = sparse(n,n);
%     IJ = ceil(n*n*rand(n,1));
%     A(IJ) = 1;
%     spy(A);
%     [V,E] = axy2ve(A,xy);
%
% Web Resources:
%   <a href="http://en.wikipedia.org/wiki/Graph_%28mathematics%29">Graphs</a>
%   <a href="http://en.wikipedia.org/wiki/Adjacency_matrix">Adjacency Matrix</a>
%
% See also: ve2axy, gplot, gplotd, gplotdc, distmat
%
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
% Release: 1.0
% Release Date: 5/22/08

error(nargchk(2,2,nargin));

[m,n] = size(A);
p = size(xy,1);
if n ~= m
    error('Invalid Adjacency Matrix: Size not square.');
elseif n ~= p
    error('Invalid XY Points: Size does not match Adjacency matrix.');
end

V = xy;
try
    [J,I] = find(A');
catch me
    try
        [I,J] = find(A);
    catch me2
        disp(me);
        rethrow(me2);
    end
end
E = [I J];
