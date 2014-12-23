function drawEigen(p, q)
  
  if (size(p,1) == size(q,1)) 
    Pts = q;
    ngray = 10;
    map = (1.0-gray(ngray))*(1.0-0.1);
    colormap(map);
    [minP  maxP] = range2(p);
    r = maxP-minP;
    if r == 0
      r = 1;
    end
    clrId = ceil(ngray * (p-minP)/1);
    clrId = max(clrId,1);
    clrId = min(clrId,ngray);
    hold on;
    for j = ngray:-1:1
      idj = (clrId == j);
      if any(idj)
        plot(Pts(idj,1), Pts(idj,2), 's','color', map(j,:),...
             'MarkerFaceColor', map(j,:),...
             'MarkerSize', 10);
      end
    end
    hold off;
    
  elseif size(p,1) == prod(q(:))
    sizeIm = q;
    showIm(reshape(p, sizeIm));
  end
  
