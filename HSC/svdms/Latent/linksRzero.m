function  [tbinNhbr] = linksRzero(A,tbinNhbr,sizeIm,R,W,figNo,flag)
%

  if ~exist('flag','var')
    flag = 1;
  end
  
  % this is how I had it but is incorrect
  %tt = find(tbinNhbr & (A > 0));
  
  % i am trying this instead 
  tt = find((A > 0));
  
  indIp = ceil(tt/prod(sizeIm)); 
  indJp = tt - (indIp-1)*prod(sizeIm); 
  [mm1,tt1] = max(W(indJp,:),[],2);
  [mm2,tt2] = max(W(indIp,:),[],2);
  RR = R(tt1 + (tt2-1)*size(R,1));
  zz = find(RR == 0);
  
  tbinNhbr(indIp(zz) + (indJp(zz)-1)*prod(sizeIm)) = 1;
  tbinNhbr(indJp(zz) + (indIp(zz)-1)*prod(sizeIm)) = 1;  
  
  
  % draw these links on the image
  indIp = indIp(zz);
  indJp = indJp(zz);
  indIx = ceil(indIp/sizeIm(1)); 
  indIy = indIp - (indIx-1)*sizeIm(1);
  indJx = ceil(indJp/sizeIm(1)); 
  indJy = indJp - (indJx-1)*sizeIm(1);
  ttx = [indIx(:) indJx(:)]';
  ttx = ttx(:);
  tty = [indIy(:) indJy(:)]';
  tty = tty(:); 
  
  if (flag & (figNo > 0))
    figure(figNo); 
    for p = 1:2:length(ttx)
      %plot(ttx(p:p+1),tty(p:p+1),'m-','linewidth',2);
      plot(ttx(p:p+1),tty(p:p+1),'b-','linewidth',3);      
    end
  end
