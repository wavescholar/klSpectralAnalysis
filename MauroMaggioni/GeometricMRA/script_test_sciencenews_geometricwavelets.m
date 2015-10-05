X=GenerateDataSets('ScienceNews');

[gW,Data,err_X,err_SVD,coeffCosts,dictCosts,lSVDidxthres] = test_threshold_mm2(X, struct('addTangentialCorrections',0,'splitting',true,'mergePsiCapIntoPhi',true,'GWTversion',1,'errorType','absolute','precision',1e-5), 1e-5:5e-5:1e-1);

figure;plot(err_X./err_SVD)

D=[cat(2,gW.WavBases{:}),cat(2,gW.ScalFuns{:})]';figure;plot(svd(D-repmat(mean(D,1),size(D,1),1)));
