function rs = infodim(s, bins)

%tstoolbox/@signal/infodim
%   Syntax:
%     * rs = infodim(s, bins)
%
%   Input arguments:
%     * s - data points (row vectors)
%     * bins - maximal number of partition per axis, default is 100
%
%   Compute the information dimension of a time-delay reconstructed
%   timeseries s for dimensions from 1 to D, where D is the dimension of
%   the input vectors. Using boxcounting approach. Scale data to be within
%   0 and 1. Give a sortiment of (integer) partitionsizes with almost
%   exponential behaviour.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

error(nargchk(1,2,nargin));

if ndim(s) ~= 2
    error('Signal must contain vector data');    
end

if nargin<2
    bins = 100;    
end

points = data(s);
[N,dim] = size(points);

% scale data to be within 0 and 1
points = points - min(min(points));
points = points / max(max(points));

% give a sortiment of (integer) partitionsizes with almost exponential behaviour 
par = [2 3 4 5 6 7 8 10 12 14 16 20 23 27 32 39 46 54 64 77 91 108 128 153 182 216 256 ...
      305 363 431 512 609 725 862 1024 1218 1449 1723 2048 2436 2897 3445 4096 4871 ...
	  5793 6889 8192 9742 11586 13778 16384 19484 23171 27555 32768 38968 46341 55109 65536];

partitions = par(find(par<=bins));	% use no sizes greater than bins

[dummy,d] = boxcount(points, partitions);

d = [zeros(1,dim) ; d];             % add zeros for partition size 1
partitions = [1 ; partitions(:)];

rs = signal(core(d), s);	
a = achse(-log2(partitions));     % create axis with arbitrary spacing
a = setname(a, 'ld r');
a2 = setname(achse(unit, 1, 1), 'Embedding dimension');
rs = setaxis(rs, 1, a);
rs = setaxis(rs, 2, a2);
rs = setplothint(rs, 'multigraph');
rs = addhistory(rs,  ['Computed information dimension']);
rs = addcommandlines(rs, 's = infodim(s', bins);
rs = setyname(rs, 'I(r)');
rs = setlabel(rs, 'Scaling of D1');
