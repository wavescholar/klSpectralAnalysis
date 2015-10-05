function  pixIds = linksKernelInter(K,rbinCut0,selectId,sizeIm,figNo,flag)

  if ~exist('flag','var')
    flag = 1;
  end
  
  
  pixIds = [];
  debug = 0;
  xids = [];
  yids = [];
  for hostId = 1:size(rbinCut0,2)
    %pixIds = kernelInter(K(:,id1),K(:,id2),sizeIm,1);
    cutIds = find(rbinCut0(:,hostId));
    if any(cutIds)
      cutIds = cutIds(cutIds > hostId);
    end
    if any(cutIds)

      for nbrId = 1:length(cutIds)
	%[selectId(hostId) selectId(cutIds(nbrId))]
	% get row/col positions of parent/nbr ids
	tt   = [selectId(hostId) selectId(cutIds(nbrId))]';
	indI = ceil(tt/sizeIm(1)); 
	indJ = tt - (indI-1)*sizeIm(1);
	xids = [xids indI'];
	yids = [yids indJ'];
	
	if (debug)
	  [ hostId cutIds(nbrId) ]
	end
	pixIdsNew = kernelInter(K(:,hostId),K(:,cutIds(nbrId)),sizeIm,debug);
	pixIds = [pixIds pixIdsNew'];
	if (debug)
	  pause;
	end
      end
      
    end
  end

  % remove duplicated pixel ids
  pixIds = unique(pixIds);
  pixIds = pixIds(:);


  % show pixIds
  if (figNo > 0)
    figure(figNo); hold on;
  end
  
  %showIm(im); hold on;
  %% their coordinates
  cols = ceil(pixIds/sizeIm(1));
  rows = pixIds - (cols - 1)*sizeIm(1);
  if (flag & (figNo > 0))
    %plot(cols,rows,'cx','linewidth',2);
    plot(cols,rows,'yx','markersize',20,'linewidth',3);
  end
  if (figNo > 0)
    title(sprintf('Pixel Count: %d',length(pixIds)));
    for p = 1:2:length(xids)
      %plot(xids(p:p+1),yids(p:p+1),'g-','linewidth',1.5);
      plot(xids(p:p+1),yids(p:p+1),'w-','linewidth',3);      
    end
  end
  
  
