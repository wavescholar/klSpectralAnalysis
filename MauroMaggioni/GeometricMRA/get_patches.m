function P = get_patches(I, patchSize, overlapping)

% Input:
% I: mxn image
% patchSize: an integer
% overlapping: a number between 0 and patchSize-1 
%                  (0=nonoverlapping, patchSize-1 = using all patches)
%
% Output: 
% P: a matrix whose rows are the image patches.

[m,n] = size(I);

spacing = patchSize-overlapping;

P = img2vecs(I, patchSize);

rows = 1:spacing:(m+1-patchSize);
cols  = 1:spacing:(n+1-patchSize);

indices = repmat(rows', 1, length(cols)) + repmat((cols-1)*(m+1-patchSize), length(rows),1);
P = P(indices,:);
