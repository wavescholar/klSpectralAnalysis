function u = yunit(d)

%tstoolbox/@description/yunit
%   return y-unit of the sampled data values (e.g. Volt, Pa etc.)
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

error(nargchk(1,1,nargin));

u = d.yunit;
