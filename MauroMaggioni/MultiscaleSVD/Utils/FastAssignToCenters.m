function vInvIdxs = FastAssignToCenters( cCenterIdxs, cPts, cDistInfo )

%
% function vInvIdxs = FastAssignToCenters( cCenterIdxs, cPts, cDistInfo )
%
% IN:
%   cCenterIdxs     : indices of points that are centers
%   cPts            : indices of points to be assigned to centers
%   cDistInfo       : DistInfo structure
%
% OUT:
%   vInvIdxs        : vInvIdxs(i) is the index into cCenterIdxs of the closest center
%                       point to cPts(i)
%

%
% SC:
%   MM  :   08/08/08
%   MM  :   08/21/09 : bug fix and speed up
%


vInvIdxs(length(cPts)) = uint32(0);

if ~iscell(cDistInfo.idxs),
    for lk = 1:length(cPts),
        for lj = 1:cDistInfo.count(lk),
            [tf,loc] = ismember( cDistInfo.idxs(lk,lj),cCenterIdxs );
            if tf,
                vInvIdxs(lk) = loc;
                break;
            end;
        end;
    end;
else
    for lk = 1:length(cPts),
        for lj = 1:cDistInfo.count(lk),
            [tf,loc] = ismember( cDistInfo.idxs{lk}(lj),cCenterIdxs );
            if tf,
                vInvIdxs(lk) = loc;
                break;
            end;
        end;
    end;
end;

return;