function [A] = shiftAffty(im,rho, minMedian)
%function [A] = shiftAffty(im,rho, minMedian)

  if nargin < 3
    minMedian = 0;
  end
%=========================
%===== Parameters  =======
%=========================
sizeIm  = size(im);
n = prod(sizeIm);
mxNhbr = 9;
  
imB = embed(im, sizeIm+2);
cx = 1+(1:sizeIm(2));
cy = 1+(1:sizeIm(1));
sxy = [1 0;  %E
       0 1;  %S
       1 1;  %SE
       -1 1];%SW
d = zeros([sizeIm 4]);
for s=1:4
  d(:,:,s) = abs(im - imB(cy+sxy(s,2), cx+sxy(s,1)));
end

%% Ignore affinities for pixels with ANY edge outside of the image when computing
%% the median.  (This is not quite all the available affinities, but it
%% is simple.
cx = 2:(sizeIm(2)-1);
cy = 1:(sizeIm(1)-1);

mdist = median(reshape(d(cy,cx,:), length(cx)*length(cy)*4, 1 ));
mdist = max(mdist, minMedian)
sigma = rho* mdist;

d = exp(-d.^2/(2*sigma^2))+eps;  %%WARNING: Need affinities >0 for binNhbr.
d = reshape(d, [prod(sizeIm) 4]);

%% Stuff d into diags of sparse affinity matrix
% Which diags?  
dgs = [sxy(:,2) + sizeIm(1)*sxy(:,1)];

% Zero out affinities for edges reaching out of image, and
% align the diags in d with the start of the matrix.
[ix iy] = meshgrid(1:sizeIm(2), 1:sizeIm(1));
ix = ix(:); iy = iy(:);
idx = iy + sizeIm(1)*(ix-1);
for s=1:4
  sx = sxy(s,1);
  sy = sxy(s,2);
  
  ids = iy+sy + sizeIm(1)*(ix+sx-1);
  ibnds = (ix+sx)<=0 | (ix+sx)>sizeIm(2) | ...
          (iy+sy)<=0 | (iy+sy)> sizeIm(1);
  %% Zero out the out of bounds elements.
  d(idx(ibnds), s) = 0;
  %% Shift the diagonal up, ready for stuffing into A
  if dgs(s) > 0

    d(:,s)  = [zeros(dgs(s),1); d(1:(end-dgs(s)), s)];
  elseif dgs(s) < 0
    strt = 1-dgs(s);
    d(:,s) = [d(strt:end, s); zeros(-dgs(s),1)];
  end
end
% Stuff it...
A = spdiags(d, dgs, prod(sizeIm), prod(sizeIm));
A = A + A';
A = spdiags(ones(prod(sizeIm),1), 0, A);

return;

%% Debug... try spdiags...
q = rand(7,1);
q'
S = spdiags([q q q q], [-2 -1 0 1], 7, 7);
S = full(S)

%% Debug... peek at some elements, they should match with 
%% a previously computed matrix B
for s = 1:4 
  [s sxy(s,:) dgs(s)]
  if (dgs(s) > 0)
    [A(1,1+dgs(s)) A(2,2+dgs(s))A(3,3+dgs(s)) A(4,4+dgs(s)) A(5,5+dgs(s))]
    d(dgs(s)+(1:5), s)'
  else
    [1, 1-dgs(s)]
    [A(1, 1-dgs(s)) A(2, 2-dgs(s))...
     A(3, 3-dgs(s)) A(4, 4-dgs(s)) ...
     A(5, 5-dgs(s))]
    [B(1, 1-dgs(s)) B(2, 2-dgs(s))...
     B(3, 3-dgs(s)) B(4, 4-dgs(s)) ...
     B(5, 5-dgs(s))]
    d((1:5), s)'
  end 
  pause;
end


  