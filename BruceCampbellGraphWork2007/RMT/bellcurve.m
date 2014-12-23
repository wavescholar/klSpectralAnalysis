%bellcurve.m
%Code 1.1 of Random Eigenvalues by Alan Edelman

%Experiment:  Generate random samples from the normal distribution.
%Observation: Histogram the random samples.
%Theory:      Falls on a bell curve.

trials=100000; dx=.2;

v=randn(1,trials);[count,x]=hist(v,[-4:dx:4]);
hold off, b=bar(x,count/(trials*dx),'y'); hold on
 
x=-4:.01:4;
plot(x,exp(-x.^2/2)/sqrt(2*pi),'LineWidth',2)
axis([-4 4 0 .45]);
