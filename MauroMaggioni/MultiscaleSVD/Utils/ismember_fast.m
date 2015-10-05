function [tf,loc] = ismember_fast(a,s)

% ismember: faster, stripped down version by MM
% Works only for one dimensional double arrays

numelA = numel(a);
numelS = numel(s);

% Initialize types and sizes.
tf = false(size(a));

% Handle empty arrays and scalars.
if numelA == 0 || numelS <= 1
    if (numelA == 0 || numelS == 0)
        return
        % Scalar A handled below.
        % Scalar S: find which elements of A are equal to S.
    elseif numelS == 1
        tf = (a == s);
        return
    end
else
    % General handling.
    [s,idx] = sort(s(:));
    if nargout <= 1
        tf = ismembc(a,s);
    else
        loc = ismembc2(a,s);
        tf = (loc > 0);
        loc(tf) = idx(loc(tf));
    end
end

return;
