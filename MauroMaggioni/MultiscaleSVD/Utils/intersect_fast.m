% intersect_fast_sorted: faster, stripped down version by MM

function [c,ia,ib] = intersect_fast(a,b,flag)

numelA = numel(a);
numelB = numel(b);
nOut = nargout;

c = [];

% Handle empty: no elements.
if (numelA == 0 || numelB == 0)
    % Predefine index outputs to be of the correct type.
    ia = [];
    ib = [];
    % Ambiguous if no way to determine whether to return a row or column.
    ambiguous = ((size(a,1)==0 && size(a,2)==0) || length(a)==1) && ...
        ((size(b,1)==0 && size(b,2)==0) || length(b)==1);
    if ~ambiguous
        c = reshape(c,0,1);
        ia = reshape(ia,0,1);
        ib = reshape(ia,0,1);
    end
elseif (numelA == 1)
    % Scalar A: pass to ISMEMBER to determine if A exists in B.
    [tf,pos] = ismember(a,b);
    if tf
        c = a;
        ib = pos;
        ia = 1;
    else
        ia = []; ib = [];
    end
elseif (numelB == 1)
    % Scalar B: pass to ISMEMBER to determine if B exists in A.
    [tf,pos] = ismember(b,a);
    if tf
        c = b;
        ia = pos;
        ib = 1;
    else
        ia = []; ib = [];
    end

    % General handling.

else
    % Switch to sort shorter list.
    if numelA < numelB
        if nOut > 1
            [a,ia] = sort(a);           % Return indices only if needed.
        else
            a = sort(a);
        end

        [tf,pos] = ismember_fast_sorted(b,a);     % TF lists matches at positions POS.

        where = zeros(size(a));       % WHERE holds matching indices
        where(pos(tf)) = find(pos);   % from set B, 0 if unmatched.
        tfs = where > 0;              % TFS is logical of WHERE.

        % Create intersection list.
        c = a(tfs);

        if nOut > 1
            % Create index vectors if requested.
            ia = ia(tfs);
            if nOut > 2
                ib = where(tfs);
            end
        end
    else
        if nOut > 1
            [b,ib] = sort(b);           % Return indices only if needed.
        else
            b = sort(b);
        end

        [tf,pos] = ismember_fast_sorted(a,b);     % TF lists matches at positions POS.

        where = zeros(size(b));       % WHERE holds matching indices
        where(pos(tf)) = find(pos);   % from set B, 0 if unmatched.
        tfs = where > 0;              % TFS is logical of WHERE.

        % Create intersection list.
        c = b(tfs);

        if nOut > 1
            % Create index vectors if requested.
            ia = where(tfs);
            if nOut > 2
                ib = ib(tfs);
            end
        end
    end
end

return;