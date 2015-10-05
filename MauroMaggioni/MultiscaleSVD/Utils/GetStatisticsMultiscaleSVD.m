function vStats = GetStatisticsMultiscaleSVD( cNets, cOpts )

%
% function vStats = GetStatisticsMultiscaleSVD( cNets, cOpts )
%
% IN:
%   cNets   : multiscale net structure constructed by FastRoughSetMultiscaleSVD.   
%   cOpts   : structur containing the following:
%               NumberOfPoints : number of data points
%               NetsOpts   : structure of options used to construct cNets
%               Lambda     : \lambda for Tychonov regularization
%               Iter       : number of iterations in Tychonov regularization
%               [J]        : computes J function statistics. Default: false.
%               [JFr]      : computes JFr function statistics. Default: false.
%               [Beta]     : compute Beta numbers statistics. Default: false.
%               [Volume]   : compute volume statistics. Default: false.
%               [ResVar]   : compute residual variance statistics. Default: false.
%
% OUT:
%   vStats  : a structure containing the following fields:
%               S                   : (number of scales)x(number of points)x(number of singular values) 
%                                       array of singular values.
%               Tychonov.S          : regularized estimate of S, as computed by HeatReg1D.
%               Tychonov.AverSlope  : average slope, for each point and each index of singular value.
%               Tychonov.DimEst     : an estimate of the intrinsic dimensionality at each point, 
%                                       computed by EstimateDim.
%               Volume              : as in PointStatisticsMultiscaleSVD
%               Opts                : options used (i.e. cOpts + defaults)
%               
%               
%
% USES:
%   PointStatisticsMultiscaleSVD, HeatReg1D, Aver_Slope, EstimateDim
%

% SC:
% (c) Copyright Mauro Maggioni
% Duke University
% 2008

if ~isfield(cOpts,'J'),           cOpts.J = false;            end;
if ~isfield(cOpts,'JFr'),         cOpts.JFr = false;          end;
if ~isfield(cOpts,'Beta'),        cOpts.Beta = false;         end;
if ~isfield(cOpts,'Volume'),      cOpts.Volume = false;       end;
if ~isfield(cOpts,'ResVar'),      cOpts.ResVar = false;       end;

% Performing point statistics (linear regression)...
vStats = PointStatisticsMultiscaleSVD( [], 1:cOpts.NumberOfPoints, cNets, struct('Lambda',40,'Iter',20,'FindLinearRegion',0,'J',cOpts.J,'JFr',cOpts.JFr,'Beta',cOpts.Beta,'Volume',cOpts.Volume,'ResVar',cOpts.ResVar) );

% Denoising singular values with Tychonov(grad^2)...
S = zeros(size(vStats.S));
Delta = cOpts.NetsOpts.Delta;
Lambda = cOpts.Lambda;
Iter = cOpts.Iter;
vStatsS = vStats.S;

for i=1:size(vStats.S,3),
    lS = squeeze(vStatsS(:,:,i))';
    %if max(max(abs(lS)))>sqrt(eps),
        S(:,:,i)=HeatReg1D(Delta, lS, Lambda,Iter)';
    %end;
end

vStats.Tychonov.S = S;
vStats.Tychonov.DimEst = [];
vStats.Opts = cOpts;

return;