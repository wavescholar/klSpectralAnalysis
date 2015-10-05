%% Go parallel with default configuration
try
    if matlabpool('size')==0,
        matlabpool
    end;
catch
end;

load MNIST_digit1

relativeprecision = [1e-1,0.5e-1,0.1e-1];
dictsizes         = [500:500:2500];

ksvdCost = zeros(length(relativeprecision),length(dictsizes));
Edata    = zeros(size(ksvdCost));
precision= zeros(size(ksvdCost));

%X0 = X0(1:550,:);

for i = 1:length(relativeprecision),
    precision = relativeprecision(i)*mean(sqrt(sum(X0.^2,2)));
    parfor k = 1:length(dictsizes),
        dictsize(i,k) = dictsizes(k);
        Edata(i,k)    = precision;
        fprintf('\n Doing (%d,%d)....',i,k);
        [D,Gamma,err] = ksvd(struct('data',X0','Edata',Edata(i,k),'dictsize',dictsize(i,k)));
        ksvdCost(i,k) = err(end)*size(X0,1)+numel(D);
        Y0 = D*Gamma;
        ksvdErr(i,k)  = norm(X0'-Y0,'fro')/size(X0,1);
        ksvdRelErr(i,k) = 0;
        for l = 1:size(X0,1),
            ksvdRelErr(i,k) = ksvdRelErr(i,k) + sqrt(sum((X0(l,:)'-Y0(:,l)).^2)/sum(X0(l,:).^2));
        end;
        ksvdRelErr(i,k) = ksvdRelErr(i,k)/size(X0,1);
        fprintf('done.');
    end;
end;




