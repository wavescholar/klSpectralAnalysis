function vDistMatrix = FastDistRestrict( cDistInfo, cIdxs )

%
% Restricts a DistInfo structure to a subset of the indices of the points
%
%

lN = length(cIdxs);

vDistMatrix = double(lN);

if ~iscell(cDistInfo.idxs),    
    for lk = 1:lN,
        lIdxs      = cDistInfo.isinidxs(cIdxs(lk),cIdxs);
        lInIdxs    = find(cDistInfo.isin(cIdxs(lk),cIdxs));
        vDistMatrix(lk,lInIdxs) = double(cDistInfo.dists(cIdxs(lk),lIdxs(lInIdxs)));
    end;
else   
    for lk = 1:lN,
        lIdxs      = cDistInfo.isinidxs{cIdxs(lk)}(cIdxs);
        [lTmp,lInIdxs]    = find(cDistInfo.isin{cIdxs(lk)}(cIdxs));
        vDistMatrix(lk,lInIdxs) = double(cDistInfo.dists{cIdxs(lk)}(lIdxs(lInIdxs)));
    end;
end;

return;



% % % if ~iscell(cDistInfo.idxs),    
% % %     for lk = 1:length(cIdxs),
% % %         lIdxs      = cDistInfo.isinidxs(cIdxs(lk),cIdxs);
% % %         lInIdxs    = find(cDistInfo.isin(cIdxs(lk),cIdxs));
% % %         vDistMatrix(lk,lInIdxs) = double(cDistInfo.dists(cIdxs(lk),lIdxs(lInIdxs)));
% % %     end;
% % % else   
% % %     for lk = 1:length(cIdxs),
% % %         lIdxs      = cDistInfo.isinidxs{cIdxs(lk)}(cIdxs);
% % %         [lTmp,lInIdxs]    = find(cDistInfo.isin{cIdxs(lk)}(cIdxs));
% % %         vDistMatrix(lk,lInIdxs) = double(cDistInfo.dists{cIdxs(lk)}(lIdxs(lInIdxs)));
% % %     end;
% % % end;
