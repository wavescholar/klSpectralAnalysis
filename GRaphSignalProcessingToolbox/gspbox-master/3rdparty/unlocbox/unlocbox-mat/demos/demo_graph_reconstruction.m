%DEMO_GRAPH_RECONSTRUCTION Reconstruction of missing sample on a graph
%
%   Here we try to reconstruct missing sample of a graph. We make the 
%   assumption that the signal is smooth. To reflect this, we chose as 
%   prior the 2 norm of the gradient of the signal. Computing the proximal 
%   operator of this function is equivalent to do a low filtering. 
%
%   For more information, please refer to references.
%
%   For this example, you need swgt toolbox. You can download it:
%   http://wiki.epfl.ch/sgwt/documents/sgwt_toolbox-1.01.zip
%
%   The problem can be expressed as this
%
%        argmin  ||Ax-b||^2 + tau*||grad(x)||_2^2
%
%   Where b is the degraded image, I the identity and A an operator
%   representing the mask. 
%
%
%
%   We set 
%
%    f_1(x)=||nabla x _2^2
%     We define the prox of f_1 as: 
%
%        prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma  ||grad(z)||_2^2
%
%     For a graph, the squared norm of the gradient is given by nabla x||_2^2 = x'Lx
%     The prox of x'Lx is:
%
%             H(x) ={1/(2*gamma*L+I}*x
%
%     We approximate this operator with a truncated Chebysheff sum. See
%     references.
%
%    f_2(x)=||Ax-b||_2^2
%     We define the gradient as: 
%
%        grad_f(x) = 2 * A^*(Ax-b)
%
%   Results
%   -------
%
%   Figure 1: Original signal on graph
%
%      This figure shows the original signal on graph.
%
%   Figure 2: Depleted signal on graph
%
%      This figure shows the signal on graph after the application of the
%      mask. Half of the vertices are set to 0.
%
%   Figure 3: Reconstructed signal on graph
%
%      This figure shows the reconstructed signal thanks to the algorithm.
%
%   References:
%     P. Combettes and J. Pesquet. A douglas-rachford splitting approach to
%     nonsmooth convex variational signal recovery. Selected Topics in Signal
%     Processing, IEEE Journal of, 1(4):564-574, 2007.
%     
%     D. Shuman, P. Vandergheynst, and P. Frossard. Chebyshev polynomial
%     approximation for distributed signal processing. In Distributed
%     Computing in Sensor Systems and Workshops (DCOSS), 2011 International
%     Conference on, pages 1-8. IEEE, 2011.
%     
%     D. Hammond, P. Vandergheynst, and R. Gribonval. Wavelets on graphs via
%     spectral graph theory. Applied and Computational Harmonic Analysis,
%     30(2):129-150, 2011.
%     
%
%   Url: http://unlocbox.sourceforge.net/doc/demos/demo_graph_reconstruction.php

% Copyright (C) 2012-2013 Nathanael Perraudin.
% This file is part of UNLOCBOX version 1.3.135
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Author: Nathanael Perraudin, Gilles Puy
% Date: sept 30 2011


%% Initialisation

clear all;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2;    % verbosity level


%% Creating the problem

% Loading data
load('minnesota.mat'); % graph
load('graph_signal.mat'); % set of value for the vertex continuous values


x=xy(:,1);
y=xy(:,2);
graph_value=s;

L=sgwt_laplacian(A); % creating the laplacian of the graph


%create the masque
p=0.5; %probability of having no label on a vertex.
M=(rand(size(graph_value))>p);

%applying the Mask to the data
depleted_graph_value=M.*graph_value;

% setting the function f2 
f2.grad=@(x) 2*M.*(M.*x-depleted_graph_value);
f2.eval=@(x) norm(M.*x-depleted_graph_value)^2;


% setting the function f1 
tau=4; 
paramsolver.gamma=0.5; % stepsize (beta is equal to 2)
h=@(x) 1./(2*paramsolver.gamma/tau.*x+1);


%Parametres for the chebysheff approximation calculaion
lmax = eigs(L,1); % approximation of the maximum eigenvalue of L
arange=[0 lmax]; % range of the eigenvalue (L is positive semi definite 
                 %=> the smallest eigenvalue is 0)
% order of the chebycheff polynomial approximation
m=10; 
N=20;

f1.prox=@(x,T) sgwt_cheby_op(x,L,sgwt_cheby_coeff(h,m,N,arange), arange);
f1.eval=@(x) tau*x'*L*x;   

%% solve the problem

% setting different parameter for the simulation
paramsolver.verbose = verbose;  % display parameter
paramsolver.maxit = 50;         % maximum number of iterations
paramsolver.tol = 10e-5;        % tolerance to stop iterating
paramsolver.method = 'ISTA';    % desired method for solving the problem

sol=forward_backward(depleted_graph_value,f1,f2,paramsolver);

%% Print the result


% Let show the original graph
msize=100;
figure(1)

set(gcf,'renderer','zbuffer');
fprintf('Displaying traffic graph\n');
set(gcf,'position',[0,600,400,400]);

hold on
scatter(x,y,msize,graph_value,'.');

[ki,kj]=find(A);
plot([x(ki)';x(kj)'],[y(ki)';y(kj)'],'k');
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);

axis equal
axis off
colorbar('location','southoutside')
caxis([0 1])
title('Original graph')


% Let show depleted graph
figure(2)
set(gcf,'renderer','zbuffer');
fprintf('Displaying the masked traffic graph\n');
set(gcf,'position',[0,600,400,400]);

hold on
scatter(x,y,msize,depleted_graph_value,'.');

[ki,kj]=find(A);
plot([x(ki)';x(kj)'],[y(ki)';y(kj)'],'k');
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);


axis equal
axis off
colorbar('location','southoutside')
caxis([0 1])
title('Deteriored graph')


% Let show the reconstructed graph
figure(3)
set(gcf,'renderer','zbuffer');
fprintf('Displaying reconstructed traffic graph\n');
set(gcf,'position',[0,600,400,400]);

hold on
scatter(x,y,msize,sol,'.');

[ki,kj]=find(A);
plot([x(ki)';x(kj)'],[y(ki)';y(kj)'],'k');
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);


axis equal
axis off
colorbar('location','southoutside')
caxis([0 1])
title('Reconstructed graph')



