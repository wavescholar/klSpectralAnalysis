function n = yname(d)

%tstoolbox/@description/yname
%   return name of the measured data (e.g. 'Heartbeat rate', 'Current'
%   etc.)
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

error(nargchk(1,1,nargin));

n = d.yname;
