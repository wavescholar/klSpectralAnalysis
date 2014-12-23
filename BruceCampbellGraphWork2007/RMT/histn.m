function [h,hn,xspan]=histn(data,x0,binsize,xf);
%HISTN Normalized Histogram.
%    [H,HN,XSPAN] = HISTN(DATA,X0,BINSIZE,XF) generates the normalized
%    histogram of area 1 from the values in DATA which are binned into
%    equally spaced containers that span the region from X0 to XF
%    with a bin width specified by BINSIZE.
%    
%     X0, BINSIZE and XF are all scalars while DATA is a vector.
%     H, HN and XSPAN are equally sized vectors.
%  
%     References:
%     [1] Alan Edelman, Handout 2: Histogramming, 
%                       Fall 2004, Course Notes 18.338.
%     [2] Alan Edelman, Random Matrix Eigenvalues.
%
%     Alan Edelman and Raj Rao, Sept. 2004.
%     $Revision: 1.1 $  $Date: 2004/09/10  17:11:18 $

xspan=[x0:binsize:xf];     

h=hist(data,xspan);          % Generate histogram
hn=h/(length(data)*binsize); % Normalize histogram to have area 1

bar(xspan,hn);               % Plot histogram
