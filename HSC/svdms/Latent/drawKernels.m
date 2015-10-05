function  drawKernels(im,selectId,sizeIm,figNo,color)

  if (nargin < 5)
    color = [1 0 0];
  end
  
  
  kCols = ceil(selectId/sizeIm(1));
  kCols = kCols(:);
  kRows = selectId(:) - (kCols - 1)*sizeIm(1);
  kText = [];
  for t = 1:length(selectId)
    kText = strvcat(kText,sprintf('%d',t));
  end
  figure(figNo); clf;
  set(figNo,'DoubleBuffer','on');  
  image(im);
  colormap(gray(256));
  resizeImageFig(figNo, size(im), 12); hold on;
  %set(get(figNo,'CurrentAxes'),'Ydir','reverse');
  axis equal;
  
  %showIm(im); hold on;
  text(kCols,kRows,kText,'color',color,'fontsize',30,'fontweight','demi');
