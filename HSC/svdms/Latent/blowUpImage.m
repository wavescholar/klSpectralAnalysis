function [sbigIm] = blowUpImage(im,scale,ranges)

  if nargin < 3
    [mn mx] = range2(im);
  else
    mn = ranges(1)
    mx = ranges(2)
  end
  %imScl = (im - mn)/(mx - mn);
  imScl = im;
  
  hb = scale;
  sizeIm = size(im);
  sbigIm = zeros(sizeIm*hb);
  for k = 1:sizeIm(1)
    for m = 1:sizeIm(2)
      sbigIm([1:hb]+(hb)*(k-1), [1:hb]+(hb)*(m-1)) = imScl(k,m);
    end
  end
