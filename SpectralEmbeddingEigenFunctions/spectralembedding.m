%% aug 4
 %playing around with spectral embedding
 
 %% load an image
 if 1 % random image
     sizebox1 = [20 20];
     sizebox2 = [8 8];
     sizebox3 = [8 8];
     mu1 = 2*ones(sizebox1); sd1 = 1;
     mu2 = 10*ones(sizebox2); sd2 = 1;
     mu3 = 6*ones(sizebox2); sd3 = 1.5;
     data1 = mu1 + sd1.*randn(sizebox1);
     data2 = mu2 + sd2.*randn(sizebox2);
     data3 = mu3 + sd3.*randn(sizebox3);
     data2border = zeros(sizebox1);
     pos2 = (sizebox1(1) - sizebox2(1))/2 + 1;
     pos3 = floor(pos2/2);
     data2border(pos3:pos3+sizebox3(1)-1,pos3:pos3+sizebox3(2)-1) = data3;
     data2border(pos2:pos2+sizebox2(1)-1,pos2:pos2+sizebox2(2)-1) = data2;
     data = data1 + data2border;
     ndata = size(data);
 else % test image
     
     data = imread('test_img_sm.jpg');
     fstart(); image(data); fstop('Nuclei Image','NuclieImage2');
     %      d = data(286:327,8:60);
     %      image(d)
     %      data = double(d); clear d;
     d = rgb2gray(data);
     d = double(d);
     d = d(2:54,:);
     data = d;
     
    
     fstart();  imagesc(d); fstop('Nuclei Image','NuclieImage');
 end
 
 % chapters 3 and 4
 % save memory:
 specdiff_mem.m % it's better to run this from the command line (open it up and run each block separately)
 
 % lots of output, does not save memory:
 specdiff_long % more output...
 
