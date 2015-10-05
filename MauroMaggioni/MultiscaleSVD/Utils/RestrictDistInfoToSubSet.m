function vDistInfo = RestrictDistInfoToSubSet( cDistInfo, cIdxs, cOpts )

%
% vDistInfo = RestrictDistInfoToSubSet( cDistInfo, cIdxs )
%
% IN:
%   cDistInfo       : a DistInfo structure
%   cIdxs           : indices of the subset to be extracted.
%   [cOpts]         : structure with the following options
%                       [RestrictOnlyRange] : restricts indices and distances of the range points, 
%                                               not of the query points. Default: false.
%                       [Idxs_sorted]       : whether cIdxs is sorted or not. Default: 0.
%                       [KnownDistanceUpperBound] : known upper bound on distances between points. Default: Inf.
%
% OUT:
%   vDistInfo       : a DistInfo structure containing only the cIdxs subset of cDistInfo. vDistInfo dists are sorted in ascending order.
%
%
% EXAMPLE:
%   lN = 100; lParam=linspace(0,2*pi,lN);cX=[sin(lParam)',cos(lParam)'];
%   [DistInfo.count,DistInfo.idxs, DistInfo.dists] = nrsearch(cX', cX', 10);
%   DistInfo2 = RestrictDistInfoToSubSet( DistInfo, 1:2:100 );
%   lOpts = struct('Delta',0.1,'DoSVD',true,'AvoidSmallClusters',0,'DistInfo',DistInfo);%,'RestrictSetIdxs',1:5:100);
%   profile on;tic;vNet = FastRoughSetDownsample( cX, lOpts );toc;profile viewer;
%   DisplayNetWithSVD( cX, vNet.Idxs, vNet.S, vNet.V );axis equal;
%   gca; hold on;plot(cX(find(vNet.InvIdxs==10),1),cX(find(vNet.InvIdxs==10),2),'k.');plot(cX(vNet.Idxs(10),1),cX(vNet.Idxs(10),2),'ro');
%

% SC:
%   MM:     4/18/2008
%   MM:     8/05/2008
%

if nargin<3,
    cOpts = [];
end;
if ~isfield(cOpts,'RestrictOnlyRange'), cOpts.RestrictOnlyRange = false; end;
if ~isfield(cOpts,'Idxs_sorted'), cOpts.Idxs_sorted = false; end;
if ~isfield(cOpts,'KnownDistanceUpperBound'), cOpts.KnownDistanceUpperBound = Inf;end;

lN = length(cIdxs);

if cOpts.KnownDistanceUpperBound<Inf,
    if iscell(cDistInfo.dists),
        vDistInfo.count = zeros(1,lN);
        vDistInfo.idxs{lN} = [];
        vDistInfo.dists{lN} = [];

        for lk = 1:length(cIdxs),
            lTmpDists = cDistInfo.dists{cIdxs(lk)};
            for lj = 1:cDistInfo.count(cIdxs(lk)),
                if lTmpDists(lj)>cOpts.KnownDistanceUpperBound,
                    break;
                end;
            end;
            lTmpDists = lTmpDists(1:lj);
            vDistInfo.idxs{lk} = zeros(1,lN);
            vDistInfo.dists{lk} = zeros(1,lN);
            [lCap,lCapIdxs,lCapIdxsIdxs] = intersect_fast_sorted_first( cIdxs, cDistInfo.idxs{cIdxs(lk)}(1:lj) );
            vDistInfo.idxs{lk} = uint32(lCapIdxs);
            vDistInfo.dists{lk} = lTmpDists(lCapIdxsIdxs);
            vDistInfo.count(lk) = length(lCap);
            [vDistInfo.dists{lk},lSortedIdxs] = sort(vDistInfo.dists{lk});
            vDistInfo.idxs{lk} = vDistInfo.idxs{lk}(lSortedIdxs);
        end;
    else
        vDistInfo.count(1,lN)  = uint32(0);
        vDistInfo.idxs(lN,lN)  = uint32(0);
        vDistInfo.dists(lN,lN) = single(0);

        for lk = 1:length(cIdxs),
            lTmpDists = cDistInfo.dists(cIdxs(lk),:);
            for lj = 1:cDistInfo.count(cIdxs(lk)),
                if lTmpDists(lj)>cOpts.KnownDistanceUpperBound,
                    break;
                end;
            end;
            lTmpDists = lTmpDists(1:lj);
            [lCap,lCapIdxs,lCapIdxsIdxs] = intersect_fast_sorted( cIdxs, cDistInfo.idxs(cIdxs(lk),1:lj) );
            vDistInfo.idxs(lk,1:length(lCapIdxs)) = uint32(lCapIdxs);
            vDistInfo.dists(lk,1:length(lCapIdxs)) = lTmpDists(lCapIdxsIdxs);
            vDistInfo.count(lk) = length(lCap);
            [vDistInfo.dists(lk,1:length(lCapIdxs)),lSortedIdxs] = sort(vDistInfo.dists(lk,1:length(lCapIdxs)));
            vDistInfo.idxs(lk,1:length(lCapIdxs)) = vDistInfo.idxs(lk,lSortedIdxs);
        end;
    end;
else
    if iscell(cDistInfo.dists),
        vDistInfo.count = zeros(1,lN);
        vDistInfo.idxs{lN} = [];
        vDistInfo.dists{lN} = [];
        if cOpts.Idxs_sorted,
            if ~cOpts.RestrictOnlyRange,
                for lk = 1:length(cIdxs),
                    vDistInfo.idxs{lk} = zeros(1,lN);
                    vDistInfo.dists{lk} = zeros(1,lN);
                    [lCap,lCapIdxs,lCapIdxsIdxs] = intersect_fast_sorted_first( cIdxs, cDistInfo.idxs{cIdxs(lk)} );
                    vDistInfo.idxs{lk} = lCapIdxs;
                    vDistInfo.dists{lk} = cDistInfo.dists{cIdxs(lk)}(lCapIdxsIdxs);
                    vDistInfo.count(lk) = length(lCap);
                    [vDistInfo.dists{lk},lSortedIdxs] = sort(vDistInfo.dists{lk});
                    vDistInfo.idxs{lk} = vDistInfo.idxs{lk}(lSortedIdxs);
                end;
            else
                for lk = 1:length(cDistInfo.idxs),
                    [lCap,lCapIdxs,lCapIdxsIdxs] = intersect_fast( uint32(cIdxs), uint32(cDistInfo.idxs{lk}) );
                    vDistInfo.idxs{lk} = uint32(lCapIdxs);
                    vDistInfo.dists{lk} = cDistInfo.dists{lk}(lCapIdxsIdxs);
                    vDistInfo.count(lk) = length(lCap);
                    [vDistInfo.dists{lk},lSortedIdxs] = sort(vDistInfo.dists{lk});
                    vDistInfo.idxs{lk} = uint32(vDistInfo.idxs{lk}(lSortedIdxs));
                end;
            end
        else
            if ~cOpts.RestrictOnlyRange,
                for lk = 1:length(cIdxs),
                    vDistInfo.idxs{lk} = zeros(1,lN);
                    vDistInfo.dists{lk} = zeros(1,lN);
                    [lCap,lCapIdxs,lCapIdxsIdxs] = intersect_fast( cIdxs, cDistInfo.idxs{cIdxs(lk)} );
                    vDistInfo.idxs{lk} = lCapIdxs;
                    vDistInfo.dists{lk} = cDistInfo.dists{cIdxs(lk)}(lCapIdxsIdxs);
                    vDistInfo.count(lk) = length(lCap);
                    [vDistInfo.dists{lk},lSortedIdxs] = sort(vDistInfo.dists{lk});
                    vDistInfo.idxs{lk} = vDistInfo.idxs{lk}(lSortedIdxs);
                end;
            else
                for lk = 1:length(cDistInfo.idxs),
                    [lCap,lCapIdxs,lCapIdxsIdxs] = intersect_fast( uint32(cIdxs), uint32(cDistInfo.idxs{lk}) );
                    vDistInfo.idxs{lk} = uint32(lCapIdxs);
                    vDistInfo.dists{lk} = cDistInfo.dists{lk}(lCapIdxsIdxs);
                    vDistInfo.count(lk) = length(lCap);
                    [vDistInfo.dists{lk},lSortedIdxs] = sort(vDistInfo.dists{lk});
                    vDistInfo.idxs{lk} = uint32(vDistInfo.idxs{lk}(lSortedIdxs));
                end;
            end
        end;
    else
        vDistInfo.count(1,lN)  = uint32(0);
        vDistInfo.idxs(lN,lN)  = uint32(0);
        vDistInfo.dists(lN,lN) = single(0);

        if cOpts.Idxs_sorted,
            if ~cOpts.RestrictOnlyRange,
                for lk = 1:length(cIdxs),
                    [lCap,lCapIdxs,lCapIdxsIdxs] = intersect( cIdxs, cDistInfo.idxs(cIdxs(lk),:) );
                    vDistInfo.idxs(lk,1:length(lCapIdxs)) = lCapIdxs;
                    vDistInfo.dists(lk,1:length(lCapIdxs)) = cDistInfo.dists(cIdxs(lk),lCapIdxsIdxs);
                    vDistInfo.count(lk) = length(lCap);
                    [vDistInfo.dists(lk,1:length(lCapIdxs)),lSortedIdxs] = sort(vDistInfo.dists(lk,1:length(lCapIdxs)));
                    vDistInfo.idxs(lk,1:length(lCapIdxs)) = vDistInfo.idxs(lk,lSortedIdxs);
                end;
            else
                for lk = 1:length(cDistInfo.idxs),
                    [lCap,lCapIdxs,lCapIdxsIdxs] = intersect_fast_sorted( cIdxs, cDistInfo.idxs(lk,:) );
                    vDistInfo.idxs (lk,1:length(lCapIdxs)) = uint32(lCapIdxs);
                    vDistInfo.dists(lk,1:length(lCapIdxs)) = cDistInfo.dists(lk,lCapIdxsIdxs);
                    vDistInfo.count(lk) = length(lCap);
                    [vDistInfo.dists(lk,1:length(lCapIdxs)),lSortedIdxs] = sort(vDistInfo.dists(lk,1:length(lCapIdxs)));
                    vDistInfo.idxs(lk,1:length(lCapIdxs)) = uint32(vDistInfo.idxs(lk,lSortedIdxs));
                end;
            end;
        else
            if ~cOpts.RestrictOnlyRange,
                for lk = 1:length(cIdxs),
                    [lCap,lCapIdxs,lCapIdxsIdxs] = intersect_fast( cIdxs, cDistInfo.idxs(cIdxs(lk),:) );
                    vDistInfo.idxs(lk,1:length(lCapIdxs)) = lCapIdxs;
                    vDistInfo.dists(lk,1:length(lCapIdxs)) = cDistInfo.dists(cIdxs(lk),lCapIdxsIdxs);
                    vDistInfo.count(lk) = length(lCap);
                    [vDistInfo.dists(lk,1:length(lCapIdxs)),lSortedIdxs] = sort(vDistInfo.dists(lk,1:length(lCapIdxs)));
                    vDistInfo.idxs(lk,1:length(lCapIdxs)) = vDistInfo.idxs(lk,lSortedIdxs);
                end;
            else
                for lk = 1:length(cDistInfo.idxs),
                    [lCap,lCapIdxs,lCapIdxsIdxs] = intersect_fast( cIdxs, cDistInfo.idxs(lk,:) );
                    vDistInfo.idxs (lk,1:length(lCapIdxs)) = uint32(lCapIdxs);
                    vDistInfo.dists(lk,1:length(lCapIdxs)) = cDistInfo.dists(lk,lCapIdxsIdxs);
                    vDistInfo.count(lk) = length(lCap);
                    [vDistInfo.dists(lk,1:length(lCapIdxs)),lSortedIdxs] = sort(vDistInfo.dists(lk,1:length(lCapIdxs)));
                    vDistInfo.idxs(lk,1:length(lCapIdxs)) = uint32(vDistInfo.idxs(lk,lSortedIdxs));
                end;
            end;
        end;
    end;
end;

return;