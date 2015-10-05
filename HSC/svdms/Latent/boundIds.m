function id = boundIds(sizeIm,w)
% return ids of all the pixels inside
% a box of size: sizeIm, ignoring all those
% pixels that are within w pixels from the
% boundary.
  
  [xx,yy] = meshgrid(1:sizeIm(2),1:sizeIm(1));
  id = find( (xx > w & xx <= sizeIm(2)-w) ...
	     .* (yy > w & yy <= sizeIm(1)-w) ...
	     );
