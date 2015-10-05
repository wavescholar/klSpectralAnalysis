function [ im ] = plot_image( x,fig )
%PLOT_IMAGE Plot image plugin for the UNLocBoX
%   Usage [ im ] = plot_image( im,fig );
%
%   Input parameters:
%         x     : Structure of data
%         fig   : Figure
%
%   Output parameters:
%         im    : Input image
%
%   This plugin display the image every iterations of an algorithm. To use
%   the plugin juste define:
%       
%       fig=figure(100);
%       param.do_sol=@(x) plot_image(x,fig);
%
%   In the structure of optional argument of the solver.
%
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/plugins/plot_image.php

% Copyright (C) 2012-2013 Nathanael Perraudin.
% This file is part of UNLOCBOX version 1.3.135
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Author: Nathanael Perraudin
% Date  : 23 septembre 2013

% select the figure
if x.iter<2
    figure(fig);
end
% display the image
imagesc(x.sol);
colormap gray;
title(['Current it: ', num2str(x.iter),'   Curr obj: ', ...
    num2str(x.curr_norm)]);
axis off         
axis image  
drawnow;
% return the image
im=x.sol;


end


