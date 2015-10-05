function DistInfo = AddSortedIdxsToDistInfo( DistInfo )

%
% This function speeds up lookup table searches in DistInfo
% It adds two fields to DistInfo, is and isinidxs
% Suppose DistInfo.idxs is a matrix (e.g. the number of nearest neighbors was fixed) - analogous definitions hold when it's a cell array.
% Fix the point k and assume the the point j is the r-th nearest neighbor of k (i.e. DistInfo.idxs(k,r)=j). Then:
%   isin(k,j) = 1;      % would be 0 if j is not amonth the nearest neighbors of k
%   isinidxs(k,j) = r;  % would be 0 if j is not amonth the nearest neighbors of k
%

lString = [];
lAll = 1:length(DistInfo.count);
if iscell(DistInfo.idxs),
    for k = length(DistInfo.count):-1:1,
        [z1,z2] = ismember(lAll,DistInfo.idxs{k}(1:DistInfo.count(k)));
        DistInfo.isin{k} = z1;
        DistInfo.isinidxs{k} = uint32(z2);
        
        if rem(k,100)==0,
            if ~isempty(lString),
                fprintf(repmat(['\b'],[1,length(lString)]));
            end;
            lString = sprintf('Sorting %d....',k);
            fprintf(lString);
        end;
    end;
else    
    for k = length(DistInfo.count):-1:1,
        [z1,z2] = ismember(lAll,DistInfo.idxs(k,1:DistInfo.count(k)));
        DistInfo.isin(k,:) = z1;
        DistInfo.isinidxs(k,:) = uint32(z2);
        if rem(k,100)==0,
            if ~isempty(lString),
                fprintf(repmat(['\b'],[1,length(lString)]));
            end;
            lString = sprintf('Sorting %d....',k);
            fprintf(lString);
        end;
    end;
end;

fprintf('\n Done.\n');

return;