function fstop(titleString,fileName)

  h = title(titleString);
            print(gcf,'-djpeg',[fileName,'.jpg']);
            close gcf; 