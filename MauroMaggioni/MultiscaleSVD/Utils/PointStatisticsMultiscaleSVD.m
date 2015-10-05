function [vPointStat,vSMean, vSStd] = PointStatisticsMultiscaleSVD( cX, cPointIdxs, cNets, cOpts )

%
% function vPointStat = PointStatisticsMultiscaleSVD( cX, cPointIdxs, cNets, cOpts )
%
% IN:
%   cX          : not needed
%   cPointsIdxs : indices of the points (into 1:M) at which to compute the statistics
%   cNets       : multiscale net, with J scales, as returned from FastRoughSetMultiscaleSVD
%   [cOpts]     : structure containing the following fields:
%                   FindLinearRegion : calls FindLinearRegion to smooth the multiscale singular values. Default: 0.
%                   SFilter          : smooths with an averaging filter of length SFilter. Default: [].
%                   Beta             : computes beta numbers of not. Default: true.
%                   ResVar           : computes residual variance or not. Default: true.
%                   Volume           : computes the volume or not. Default: true;
%                   JFr              : computes the Frobenius Jones' function. Default: true.
%                   J                : computes the Jones' function. Default: true.
%
%
% OUT:
%   vPointStat  : structure containing
%                   J           : the Jones' function evaluated at the points of index cPointsIdxs
%                   JFr         : the Frobenius Jones' function evaluated at the points of index cPointsIdxs
%                   Beta        : J by M by N tensor containing the (j,x,n) Jones' beta number
%                   S           : J by M by N tensor containing the (j,x,n) singular value
%                   ResVar      : J by M by N tensor containing the sum of squares of the (j,x,n:N) singular values squared.
%                   Volume      : J by M matrix containing the volume of the Voronoi cell at scale j around each point
%                   S_sm        : denoised (across scales) S (if cOpts.SFilter or cOpts.FindLinearRegion)
%                   B_sm        : regression coefficients (if cOpts.FindLinearRegion)
%                   S_sm_Indx   : cell array of same size as S_sm with indices of the points not used in linear regression
%
%
% USES:
%   PointIndicesInMultiscaleNet
%

% SC:
%   MM:     3/10/2008   : revised
%   MM:     8/04/2008   : revised
%
% Mauro Maggioni
% (c) Copyright Duke University, 2006
%

if nargin<4,
    cOpts = [];
end;
if ~isfield(cOpts,'SFilter'),               cOpts.SFilter = [];             end;
if ~isfield(cOpts,'FindLinearRegion'),      cOpts.FindLinearRegion = 0;     end;
if ~isfield(cOpts,'Beta'),                  cOpts.Beta = true;              end;
if ~isfield(cOpts,'ResVar'),                cOpts.ResVar = true;            end;
if ~isfield(cOpts,'J'),                     cOpts.J = true;                 end;
if ~isfield(cOpts,'JFr'),                   cOpts.JFr = true;               end;


vPointStat = []; vSMean = []; vSStd = [];

lPointNetIdxs = PointIndicesInMultiscaleNet( cPointIdxs, cNets );

lDim = size(cNets(1).S,2);
lN = length(cPointIdxs);

% Compute multiscale statistics
if cOpts.J,
    vPointStat.J     (lN,lDim)               = single(0);
    vPointStat.J_lbl                         = 'J fcn';
end;
if cOpts.JFr,
    vPointStat.JFr   (lN,lDim)               = single(0);
    vPointStat.JFr_lbl                       = 'J fcn res';
end;
vPointStat.S     (length(cNets),lN,lDim) = single(0);
vPointStat.S_lbl                         = 'Sing. vals.';
if cOpts.Beta,
    vPointStat.Beta     (length(cNets),lN,lDim) = single(0);
    vPointStat.Beta_lbl                      = 'Beta';
end;
if cOpts.ResVar,
    vPointStat.ResVar(length(cNets),lN,lDim) = single(0);
    vPointStat.ResVar_lbl                    = 'Res. var.';
end;
if cOpts.Volume,
    vPointStat.Volume(length(cNets),lN)      = uint32(0);
    vPointStat.Volume_lbl                    = 'Vol.';
end;

for lj = 1:length(cNets),
    if cOpts.Volume,
        vPointStat.Volume(lj,:)          = uint32(cNets(lj).count(lPointNetIdxs(:,lj)));
    end;
    vPointStat.S(lj,:,:)                 = single(cNets(lj).S(lPointNetIdxs(:,lj),:));
    for ld = 1:lDim,
        if cOpts.J,
            vPointStat.J(:,ld)           = vPointStat.J(:,ld) + cNets(lj).S(lPointNetIdxs(:,lj),ld).^2;
        end;
        if cOpts.JFr,
            vPointStat.JFr(:,ld)         = vPointStat.JFr(:,ld) + sum(cNets(lj).S(lPointNetIdxs(:,lj),ld:lDim).^2,2);
        end;
        if cOpts.ResVar,
            vPointStat.ResVar(lj,:,ld)   = single(sum(vPointStat.S(lj,:,ld:lDim).^2,3));
        end;
        if cOpts.Beta,
            %vPointStat.Beta(lj,:,ld)     = squeeze(vPointStat.ResVar(lj,:,ld))./(double(vPointStat.Volume(lj,:))*((cNets(lj).Delta).^2));
            vPointStat.Beta(lj,:,ld)     = single(squeeze(vPointStat.ResVar(lj,:,ld))./(((cNets(lj).Delta).^2)));
        end;
    end;
end;

if ~isempty(cOpts.SFilter),
    vPointStat.S_sm(size(vPointStat.S,1),size(vPointStat.S,2)) = single(0);
    %vPointStat.B_sm(size(cX,2),size(vPointStat.S_sm,2),size(vPointStat.S_sm,3)) = single(0);
    vPointStat.B_sm = 0;
    %    fprintf('000000');
    for lk = 1:size(vPointStat.S,2),
        %        fprintf('\b\b\b\b\b\b%.6d',lk);
        for ld = 1:size(vPointStat.S,3),
            lTmp = conv((squeeze(vPointStat.S(:,lk,ld))),cOpts.SFilter);
            vPointStat.S_sm(:,lk,ld) = lTmp((round(length(cOpts.SFilter)/2)):length(lTmp)-round(length(cOpts.SFilter)/2)+1);
        end;
    end;
elseif cOpts.FindLinearRegion,
    vPointStat.S_sm = zeros(size(vPointStat.S));
    lDeltas = cat(1,cNets.Delta);
    %    fprintf('000000');
    for lk = 1:size(vPointStat.S,2),
        %        fprintf('\b\b\b\b\b\b%.6d',lk);
        for ld = 1:size(vPointStat.S,3),
            [vPointStat.S_sm(:,lk,ld),vPointStat.B_sm(:,lk,ld),vPointStat.S_sm_Indx{lk,ld}] = FindLinearRegion( lDeltas,squeeze(vPointStat.S(:,lk,ld)) );
        end;
    end;
else
    vPointStat.S_sm = [];
end;


return;