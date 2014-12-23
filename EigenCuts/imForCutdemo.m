function [im,imName] = imForCutdemo()

fNameHdr(1,:) = 'random.206863';
fNameHdr(2,:) = 'random.213440';
fNameHdr(3,:) = 'random.208346';
fNameHdr(4,:) = 'random.208212';
fNameHdr(5,:) = 'random.207747';
fNameHdr(6,:) = 'random.211072';
fNameHdr(7,:) = 'random.212814';
fNameHdr(8,:) = 'random.209857';
fNameHdr(9,:) = 'random.206152';
fNameHdr(10,:) = 'random.206153';
t = -1;
while (t < 0 |  t > 12)
    clc;
    fprintf('\nSelect one of the following:\n');
    fprintf('\t  0: create a random image \n');
    for k = 1:9
        fprintf('\t  %d: %s.pgm \n',k,fNameHdr(k,:));
    end
    fprintf('\t 10: eye.pgm \n');
    fprintf('\t 11: id25_reye016.pgm \n');
    fprintf('\t 12: test_img_sm.pgm \n');
    
    t = input('Selection> ');
    if (t < 0 | t > 11) clc; end
end

if (t)
    
    if (t <= 9)
        fName = [fNameHdr(t,:) '.pgm'];
        im = pgmRead(fName);
        imName = fNameHdr(t,:);
        
        figure(201); showIm(im);
        fprintf('Idea is to separate the square in the center ');
        fprintf('from the background..\n');
        fprintf('Hit enter to continue\n');
        pause;
        clc;
        
    else
        if (t == 10)
            im = pgmRead('eye.pgm');
            imName = 'eye';
        end
        if (t == 11)
            im = pgmRead('id25_reye016.pgm');
            imName = 'id25_reye016';
        end
        if (t == 12)
            im = pgmRead('test_img_sm.pgm');
            imName = 'test_img_sm';
        end
        figure(201); showIm(im);
        fprintf('Idea is to partition the eye image such that,\n');
        fprintf('there is a small group of coherent regions. \n');
        fprintf('Hit enter to continue\n');
        pause;
        clc;
        
        return;
    end
    
else
    
    ok = 1;
    ctr = 1;
    while (ok)
        
        seed = round(sum(100*clock));
        rand('seed', seed);
        
        fnameRoot = 'random';
        sigmar    = 3;
        
        %%% Build Gaussian filter masks, along with derivatives.
        sigmaSqr  = sigmar*sigmar;
        gFiltSize = 2 * round(3.0 * sigmar) + 1;
        x = [1:gFiltSize] - round((gFiltSize+1)/2);
        gFilt   = exp(- x .* x / (2.0*sigmaSqr));
        gFilt   = gFilt/ sum(gFilt(:));
        
        sz   = 32;
        %imo  = rand(sz,sz);
        sigman = 2.5;
        %sigman = 1.5;
        imo  = 10 + sigman*randn(sz,sz);
        imor = rconv2sep(imo, gFilt, gFilt);
        imor = imor(1:2:end, 1:2:end);
        %imor= corrDn(imo, gFilt, 'reflect1');
        %imor= corrDn(imor, gFilt', 'reflect1',[2 2]);
        
        %imi  = rand(sz,sz);
        imi  = 9.5 + sigman*randn(sz,sz);
        imir = rconv2sep(imi, gFilt, gFilt);
        imir = imir(1:2:end,1:2:end);
        
        im = imor;
        
        lo = 7;
        hi = 11;
        im(lo:hi,lo:hi) = imir(lo:hi,lo:hi);
        
        
        pgmWrite(im,'data/random.tmp.pgm');
        im = pgmRead('data/random.tmp.pgm');
        imName = [fnameRoot '.' sprintf('%d',seed) '.pgm'];
        
        figure(201); showIm(im);
        if (ctr == 1)
            clc;
            fprintf('Idea is to separate the square in the center from the\n');
            fprintf('background. The generative process satisfies the \n');
            fprintf('following criterion:\n');
            fprintf(' - provide variability in the affinities within \n');
            fprintf('   background and foreground regions\n');
            fprintf(' - there be a significant step in gray-level\n');
            fprintf('   across the boundary between the two regions\n');
            fprintf(' - there be a significant chance for pixels in \n');
            fprintf('   the opposite sides of the boundary to have\n');
            fprintf('   similar gray-levels and thus high affinities.\n');
            fprintf('\nHowever if the square is barely visible \n');
            fprintf('it may be a good idea to generate a new image \n');
        end
        
        ok = input('Enter a non-zero number for a new image> ');
        ctr = ctr + 1;
    end
    
end
  
