n=50000;
trials=100;
delta=.05;
x=-4:delta:4;
bell=exp(-x.^2/2)/sqrt(2*pi);
cov=zeros(length(x));

hold off
for i=1:trials,
    a=randn(1,n);
    y=hist(a,x);
    y=(y/n)/delta;
    z=sqrt(n)*(y-bell)./sqrt(bell);
    %plot(x, z);
    cov=cov+z'*z;
     
end
cov=delta*cov/trials;
cov=cov(2:end-1,2:end-1);
hold off
plot(diag(cov));
hold on
plot(diag(cov,1),'r');