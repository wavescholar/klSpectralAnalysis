clear all;
[X, GWTopts, imgOpts] = GenerateData_and_SetParameters('Cosine');
GWT = GMRA(X, GWTopts);
[GWT, Data] = GWT_trainingData(GWT, X);
GWT_displayResults(GWT,Data,imgOpts);
