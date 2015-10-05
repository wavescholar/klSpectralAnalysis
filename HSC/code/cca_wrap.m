function [Ad indices] = cca_wrap(A,tau,takeall)

if nargin < 2
    tau = 0;
end
if nargin < 3
    takeall = 1; %otherwise only report clusters that are at least 33% as large as maximum
end

% get largest connected component from A
% addpath('/home/chakra/virginia/Desktop/cryoEM/MEX-cca-yalla');

lbl = gingcca_symmetric(sparse(A),tau); %[lbl z] = cca(sparse(A),tau,1); (mex version is faster, but sometimes causes matlab to crash)
ulbl = unique(lbl);
nclasses = size(ulbl,1) % total number of connected components

% sort the label information
ulblLen = zeros(length(ulbl),1);
for t = 1:length(ulbl)
  ulblLen(t) = length(find(lbl == ulbl(t)));
end
[sulblLen,iulblLen] = sort(-ulblLen);
sulblLen = -sulblLen;
% look at largest indices
sulblLen(1:min(nclasses,12))'
maxSize = sulblLen(1)
if ~takeall % only return large clusters
    dividefurther_indices = iulblLen(sulblLen >= .3333*maxSize);
else
    dividefurther_indices = iulblLen;
end

ndiv = size(dividefurther_indices,1)
% does the -1 belong here??? why??!!
disp(' why is this necessary??!!');
subtract1 = 0;
% make sure the indicing is correct
dfA = A(find(lbl == dividefurther_indices(1)));
if size(dfA,1) ~= sulblLen(1)
    if ~takeall
        dividefurther_indices = iulblLen(sulblLen >= .3333*maxSize) - 1;
    else
        dividefurther_indices = iulblLen - 1;
    end
    ndiv = size(dividefurther_indices,1);
    subtract1 = 1;
    disp('subtracted 1 from dividefurther_indices');
end


ndiv = size(dividefurther_indices,1)
indices = cell(ndiv,1);
nind = zeros(ndiv,1);
for k = 1:ndiv
indices{k} = find(lbl == dividefurther_indices(k));
nind(k) = size(indices{k},1);
end
[a b] = max(nind)
Ad = A(indices{b},indices{b});
