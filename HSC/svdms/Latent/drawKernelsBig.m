function  drawKernelsBig(im,selectId,sizeIm,figNo,color)

  if (nargin < 5)
    color = [1 0 0];
  end
  
  bigIm = blowUpImage(im,8);
  
  
  kCols = ceil(selectId/sizeIm(1));
  kCols = kCols(:);
  kRows = selectId(:) - (kCols - 1)*sizeIm(1);
  kText = [];
  for t = 1:length(selectId)
    kText = strvcat(kText,sprintf('%d',t));
  end
  
  kCols = kCols*8;
  kRows = kRows*8;
  
  figure(figNo); clf;
  image(bigIm);
  resizeImageFig(figNo, size(bigIm), 2);
  %showIm(bigIm); hold on;
  text(kCols,kRows,kText,'color',color,'fontsize',12,'fontweight','demi');

