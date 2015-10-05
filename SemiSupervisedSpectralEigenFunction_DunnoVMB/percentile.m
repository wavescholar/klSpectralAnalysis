function [PerLower,PerUpper] = percentile(X,LOWER,UPPER);

  %% function that returns thresholds for LOWER and UPPER fractions of
  %the data (not percentiles, since LOWER and UPPER are 0...1 rather than 0..100).
  
  if nargin==1
    LOWER = 0.025;
    UPPER = 0.975;
  end
  
  
Y = ones(1,length(X));

%Sort them
[X, Index] = sort(X);
Y = Y(Index);

%Find the desired indices
NumSamples = sum(Y);
MedianInd = (NumSamples+1)/2;
PerLowerInd = LOWER * (NumSamples + 1);
PerUpperInd = UPPER * (NumSamples + 1);

%cumulative frequency
CumFreq = cumsum(Y);

%Get the frequency bins of interest
MedBinFloor = find(CumFreq >= floor(MedianInd), 1);
MedBinCeil = find(CumFreq >= ceil(MedianInd), 1);

PerLowerBinFloor = find(CumFreq >= floor(PerLowerInd), 1);
PerLowerBinCeil = find(CumFreq >= ceil(PerLowerInd), 1);

PerUpperBinFloor = find(CumFreq >= floor(PerUpperInd), 1);
PerUpperBinCeil = find(CumFreq >= ceil(PerUpperInd), 1);

Median = (X(MedBinFloor) + X(MedBinCeil))/2;

%Average is used instead of interpolation to match PRCTILE
PerLower = (X(PerLowerBinFloor) + X(PerLowerBinCeil))/2;
PerUpper = (X(PerUpperBinFloor) + X(PerUpperBinCeil))/2;

