function f_efunc = funcseg(xx,is,dd2,uu2) 
% adapted from fergus' demo: more pictures, more general, binning slightly
% more stable (Wsquig), also sigma not constant, except for exact solution
 
nPoints = prod(size(xx));
ndata = size(xx);
M = size(is,2); % number of labels..


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Plot data and labeled points
figure; imagesc(xx);
Colors = rand(M,3);
hold on;
for k = 1:M
    [a b] = ind2sub(ndata,is(:,k));
    plot(b,a,'mo','MarkerSize',10,'MarkerFaceColor',Colors(k,:));
end
title('are the labels okay?','FontSize',20);
pause
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% setup weights on datapoints (lambda) and labels (y)
lambda=zeros(nPoints,1);
y = zeros(nPoints,1);
%labels = [+1 -1]; % need to change if more than 2 classes!
labels = 2*[1:M]-1
for j = 1:M % for each class
    for k = 1:size(is,1)
        if is(k,j) ~= 0 % sometimes is contains 0s, if one class does not have as many labeled points as other classes
            y(is(k,j))=labels(j);
            lambda(is(k,j))=1000;
        end
    end
end

 % build diagonal Lambda matrix
 Lambda=diag(lambda);
 

 
 
 %%% now solve for coefficients in NUM_EVECS x NUM_EVECS linear system
 %%% this is eqn. 1 in the NIPS paper (but with approx. eigenvectors/values).
 alpha2=(dd2 +uu2'*Lambda*uu2)\(uu2'*Lambda*y);
 f_efunc=uu2*alpha2;
 
% multiple labels:
 cls = cell(M,1);
 cls{1} = find(f_efunc < 2);
 for k = 2:M-1
     b = f_efunc >= 2*(k-1);
     c = f_efunc < 2*k;
     cls{k} = find(b & c) ;
 end
 cls{M} = find(f_efunc >= 2*(M-1));
 figure; subplot(2,1,2); plot(f_efunc); title('labeling','Fontsize',16);
 
 subplot(2,1,1);
 xx_class = zeros(size(xx));
 for k = 1:M
    xx_class(cls{k}) = k;
 end
 imagesc(xx_class); colormap('cool')
 axis equal;
 axis off;
 title('Eigenfunction','Fontsize',16);
