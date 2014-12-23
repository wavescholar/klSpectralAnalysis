%classicexpt.m
%Code 1.2 of Random Eigenvalues by Alan Edelman and Raj Rao

%Experiment:  Generate many realizations of the classical random matrices. 
%Observation: Histogram the eigenvalues of the random matrices.
%Theory:      Falls on the predicted curves.
%Reference:   [1] Alan Edelman, Handout 3, Experiments with Classical
%                               Matrix Ensembles, Fall 2004.
%             [2] Alan Edelman, Random Matrix Eigenvalues.

trials = 1000;  % Number of trials
e = [];         % Array for collecting eigenvalues
syms f x;       % Symbolic variables for theoretical predictions
                % Requires Symbolic Toolbox in MATLAB

figure; set(gcf,'DoubleBuffer','On'); 
warning off;
for i = 1 : trials 
    
    % Generate random matrix (uncomment the appropriate line)
    
    isreal = 1;
    % N=50; C = wigner(N,isreal); f = wigneredf;
    % N=50; M=100; C = wishart(N,M,isreal); f = wishartedf(N/M);
    % N=50; M1=100; M2=100; C = manova(N,M1,M2,isreal); f = manovaedf(N/M1,N/M2);
    
    % Collect the eigenvalues 
    
    e = [e;real(eig(C))]; % Take real part to compensate for numerical issues
    
    % Plot density of eigenvalues
    
    x0=-2; binsize=0.05; xf=4;
    [h,hn,xspan] = histn(e,x0,binsize,xf); 
    axis([x0 xf 0 1.5]); hold on

    % Plot theory
     
    F = real(subs(f,xspan)); F(F==inf)=0;
    hold on
    ts=['# trials = ' num2str(i)]; title(ts,'FontSize',15);
    plot(xspan,F,'red','LineWidth',2);drawnow;
    hold off

end
