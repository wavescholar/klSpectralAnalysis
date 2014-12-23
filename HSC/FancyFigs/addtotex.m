%addtotex add the code in the tex file for each figure
%
%   fid         file id
%   figpath     string used in the argument of "includegraphics"
%   fname       string for the label
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       Jonathan Rossel     22.07.08

function addtotex(fid,figpath,fname);

%append the code
fprintf(fid,'%%%s\n',fname);
fprintf(fid,'\\begin{figure}[!htb]\n');
fprintf(fid,'\\begin{center}\n');
fprintf(fid,'\t\\includegraphics[width=\\textwidth]{%s}\n',figpath);
fprintf(fid,'\t\\captionsm{captiontext}\n');
fprintf(fid,'\t\\label{fig:%s}\n',fname);
fprintf(fid,'\\end{center}\n');
fprintf(fid,'\\end{figure}\n\n');