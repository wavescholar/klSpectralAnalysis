function em = embed(im, sz)
% em = embed(im, sz)
% embed a smaller image centered in a larger image of size sz and
% pad the surround with zeros.
 em = zeros(sz);
 im_sz = size(im);
 if (min(sz - im_sz) <0)
  error('EMBED: image larger than embedding size');
 end
 shift = round((sz - im_sz)/2);
 em((shift(1)+1):(shift(1)+im_sz(1)), (shift(2)+1):(shift(2)+im_sz(2))) = im;
 return;   
