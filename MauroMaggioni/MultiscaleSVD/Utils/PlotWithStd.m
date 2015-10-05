function [vPlotHandle,vMeans,vStds] = PlotWithStd( cX, cY, cPlotOpts )

%
% function vPlotHandle = PlotWithStd( cX, cY, cPlotOpts )
%
% IN:
%   cX      : N vector parametrizing the X axis
%   cY      : M by N matrix of M experiments and N trials for each experiment, or
%             M by N by P matrix of M experiments and N trials for P different variables
%   [cPlotOpts] : as the options for the Matlab function errorbar. Default: 'k-'.
%
% OUT:
%   vPlotHandle : handle to the object plotted
%
% SC:
%   MM  : 12/27/2005
%
% (c) Copyright 2005
%
% Mauro Maggioni (mauro.maggioni@yale.edu)
%
% Department of Mathematics
% Box 208283
% Yale University, CT 06520
%
%

if nargin<3,
    cPlotOpts = [];
end;

if (nargin<2) | (isempty(cPlotOpts)),
    cPlotOpts = 'k-';
end;

if length(size(cY)==2),
    % Compute mean and stardard deviation
    lMean = mean(double(cY),2);
    lStd  = std(double(cY),0,2);
    
    %vPlotHandle = errorbar(cX,lMean,lStd,cPlotOpts);
    
    vPlotHandle = [];
    vPlotHandle = plot(cX,lMean+lStd,'b--','LineWidth',0.5,'MarkerSize',0.5); hold on;
    vPlotHandle = [vPlotHandle;plot(cX,lMean-lStd,'b--','LineWidth',0.5,'MarkerSize',0.5)];
    vPlotHandle = [vPlotHandle;plot(cX,lMean,'k','LineWidth',1)];
    
    vMeans = lMean;
    vStds  = lStd;
else
    for p = 1:size(cY,3),
        % Compute mean and stardard deviation
        lMean = mean(double(squeeze(cY(:,:,p))),2);
        lStd  = std(double(squeeze(cY(:,:,p))),0,2);
        
        vPlotHandle = errorbar(cX,lMean,lStd,cPlotOpts); hold on;
        
        vMeans(:,p) = lMean;
        vStds(:,p)  = lStd;
    end;
end;

return;