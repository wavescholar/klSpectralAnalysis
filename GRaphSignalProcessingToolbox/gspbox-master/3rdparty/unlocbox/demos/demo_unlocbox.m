%DEMO_UNLOCBOX  Simple tutorial for the UNLocBoX
%
%   Introduction
%   ------------
%
%   This toolbox is designed to solve convex optimization problems of the
%   form:
%
%               argmin_x  f1(x) + f2(x),
%
%   or more generally
%   
%               argmin_x sum_{n=1}^K f_n(x), 
%
%   where the f_i are  lower semi-continuous convex functions and x the
%   optimization variables. For more details about the problems, please
%   refer to the userguide (UNLocBoX-note-002) availlable on
%   http://unlocbox.sourceforge.net/note/. 
%
%   This toolbox uses proximal splitting methods. Those cut the
%   problem into smaller (and easier) subproblems that will be solved
%   iteratively.
%
%   The toolbox is essentially made of three kind of functions:
%   
%    Solvers: the core of the toolbox
%
%    Proximity operators: they solve small minimization problems and 
%     allow a quick implementation of common problems. 
%
%    Demonstration files: examples to help you to use the toolbox
%
%   This toolbox is provided for free. We would be happy to receive
%   comments, bugs or any kind of help that could help us to improve the
%   toolbox.
%
%   The problem
%   -----------
%   Let's suppose we have a noisy image with missing pixels. Our goal would
%   be to find the closest image to the original one. We begin first by
%   setting up some assumptions about the problem.
%
%   Figure 1: This figure shows the original image (The cameraman).
%
%      
%
%   Assumptions
%   -----------
%   In this particular example, we firstly assume that we know the position 
%   of the missing pixels. This can be the result of a previous process on 
%   the image or a simple assumption.
%   Secondly, we assume that the image is not special. Thus, it is
%   composed of well delimited patches of colors.
%   Thirdly, we suppose that known pixels are subject to some Gaussian
%   noise with a variance of epsilon.
%
%   Figure 2: This figure shows the noisy image.
%
%      
%
%   Figure 3: This figure shows the measurements. 50 percents of the pixels have been removed.
%
%       
%
%   Formulation of the problem
%   --------------------------
%   At this point, the problem can be expressed in a mathematical form. We
%   will simulate the masking operation by a mask A. This first
%   assumption leads to a constraint.  
%
%              Ax = y
%
%   where x is the vectorized image we want to recover, y are the
%   observe noisy pixels and A a linear operator representing the mask.
%   However due to the addition of noise this constraint is a little bit
%   relaxed. We rewrite it in the following form 
%
%      .       || Ax - y ||_2  <=  epsilon
%
%   Note that epsilon can be chosen equal to 0 to satisfy exactly the 
%   constraint.
%   In our case, as the measurements are noisy, we set epsilon to the
%   standard deviation of the noise.
%
%   Using the prior assumption that the image has a small TV-norm (image
%   composed of patch of color and few degradee), we will express the
%   problem as
%
%       argmin ||x||_TV s.t ||Ax-y||_2 < epsilon        (Problem I)
%
%   where b is the degraded image and A an linear operator
%   representing the mask. epsilon is a free parameter that tunes the
%   confidence to the measurements.
%   This is not the only way to define the problem. We could also write:
%
%       argmin ||Ax-y||_2 + lambda  ||x||_TV            (Problem II)
%
%   with the first function playing the role of a data fidelity term and
%   the second a prior assumption on the signal. lambda adjusts the tradeoff
%   between measurement fidelity and prior assumption. We call it the 
%   regularization parameter. The smaller it is, the more we trust the
%   measurements and vice-versa. epsilon play a similar role as lambda.
%
%   Note that there exist a bijection between the parameters lambda and
%   epsilon leading to the same solution. The bijection function is not
%   trivial to determine. Choosing between one or the other problem will
%   affect the solvers and the convergence rate.
%
%   Proximity operator
%   ------------------
%
%   The UNLocBoX is using proximal splitting techniques for solving convex
%   optimization problem. Those techniques consist in dividing the problem
%   into easier problems. Each function is minimized with its proximity
%   operator.
%   In this particular case, the solver will iteratively minimize the TV 
%   norm and then perform the projection. 
%
%   The proximity operator of a lower semi-continuous convex function f is
%   defined by:
%
%        prox_{f,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma f(x)
%
%   Proximity operators are minimizing a function without going too far
%   from a initial point. They can be assimilated as denoising operator. They
%   are also considered as generalization of projection.
%
%   Indeed, the proximity operator of the indicator function define by:
%
%                  /   0       if   x in C   
%         i_C(x) = |
%                  \  inf      otherwise
%
%   is simply the projection onto the set C.
%
%   The toolbox solvers cannot directly handle constraints and only take
%   functions as input. However the constraints can be inserted in the 
%   objective function in the form of indicator functions. 
%
%   Solving problem I
%   -----------------
%
%   The UNLocBoX solvers take as input functions with their proximity
%   operator or with their gradient. In the toolbox, functions are
%   modelized with structure object with at least two fields. One field
%   contains an operator to evaluate the function and the other allows to
%   compute either the gradient (in case of differentiable function) or the
%   proxitity operator ( in case of non differentiable functions). In this
%   example, we need to provide two functions: 
%
%    f_1(x)=||x||_{TV}
%     We define the prox of f_1 as: 
%
%        prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma ||z||_TV
%
%     This function is defined in Matlab using:
%
%           paramtv.verbose=1;
%           paramtv.maxit=50;
%           f1.prox=@(x, T) prox_tv(x, T, paramtv);
%           f1.eval=@(x) tv_norm(x);   
%
%     This function is a structure with two fields. First, f1.prox is an
%     operator taking as input x and T and evaluating the proximity
%     operator of the function (T plays the role of gamma is the
%     equation above). Second, and sometime optional, f1.eval is also an
%     operator evaluating the function at x.
%
%     The proximal operator of the TV norm is already implemented in the
%     UNLocBoX by the function prox_tv. We tune it by setting the maximum
%     number of iterations and a verbosity level. Other parameters are also
%     available (see documentation).
%
%      paramtv.verbose selects the display level (0 no log, 1 summary at
%       convergence and 2 display all steps).
%
%      paramtv.maxit defines the maximum number of iteration.
%
%    f_2 is the indicator function of the set S defines by Ax-b||_2 <epsilon 
%     We define the proximity operator of f_2 as 
%
%        prox_{f2,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma i_S( x ),
%
%     with i_S(x) is zero if x is in the set S and infinite otherwise.
%     This previous problem has an identical solution as:
%
%        argmin_{z} ||x - z||_2^2   s.t.  || A z - y||_2 < epsilon
%
%     It is simply a projection on the B2-ball. In matlab, we write:
%
%           param_proj.epsilon=epsilon;
%           param_proj.A=A;
%           param_proj.At=A;
%           param_proj.y=y;
%           f2.prox=@(x,T) proj_b2(x,T,param_proj);
%           f2.eval=@(x) eps;
%
%     The prox field of f2 is in that case the operator
%     computing the projection. Since we suppose that the constraint is
%     satisfied, the value of the indicator function is 0. For
%     implementation reasons, it is better to set the value of the operator
%     f2.eval to eps than to 0. Note that this hypothesis could lead
%     to strange evolution of the objective function. Here the parameter
%     A and At are mendatory. Note that A = At, since the masking
%     operator can be performed by a diagonal matrix containing 1's for
%     observed pixels and 0's for hidden pixels.
%
%   At this point, a solver needs to be selected. The UNLocBoX contains many
%   different solvers. You can try them and observe the convergence speed. 
%   Just remember that some solvers are optimized for specific problems. 
%   In this tutorial, we present two of them forward_backward and
%   douglas_rachford. Both of them take as input two functions (they have
%   generalization taking more functions), a starting point and some
%   optional parameters. 
%
%   In our problem, both functions are not smooth on all points of the
%   domain leading to the impossibility to compute the gradient. In that
%   case, solvers (such as forward_backward) using gradient descent
%   cannot be used. As a consequence, we will use douglas_rachford
%   instead. In matlab, we write:
%
%           param.verbose=1;    
%           param.maxit=100;    
%           param.tol=10e-5;    
%           param.gamma=1;     
%           sol = douglas_rachford(y,f1,f2,param); 
%
%    param.verbose selects the display level (0 no log, 1 summary at
%       convergence and 2 display all steps).
%
%    param.maxit defines the maximum number of iteration.
%
%    param.tol is stopping criterion for the loop. The algorithm stops if
%
%         (  n(t) - n(t-1) )  / n(t) < tol,
%      
%     where  n(t) is the objective function at iteration t*
%
%    param.gamma defines the stepsize. It is a compromise between
%     convergence speed and precision. Note that if gamma is too big, the
%     algorithm might not converge.
%
%
%   Figure 4: This figure shows the reconstructed image by solving problem I using Douglas Rachford algorithm.
%
%          
%
%   Solving problem II
%   ------------------
%
%   Solving problem II instead of problem I can be done with a small
%   modification of the previous code. First we define another function as
%   follow:
%
%           param_l2.A=A;
%           param_l2.At=A;
%           param_l2.y=y;
%           param_l2.verbose=1;
%           f3.prox=@(x,T) prox_l2(x,lambda*T,param_l2);
%           f3.grad=@(x) 2*lambda*A(A(x)-y);
%           f3.eval=@(x) lambda*norm(A(x)-y,'fro');
%
%   The structure of f3 contains a field f3.grad. In fact, the l2 norm 
%   is a smooth function. As a consequence the gradient is well defined on 
%   the entire domain. This allows using the forward_backward solver.
%   However, we can in this case also use the douglas_rachford solver.
%   For this we have defined the field f3.prox.
%
%   We remind that forward_backward will not use the field f3.prox and
%   douglas_rachford will not use the field f3.grad.
%
%   The solvers can be called by:
%   
%           sol21 = forward_backward(y,f1,f3,param);
%
%   Or:   
%   
%           sol22 = douglas_rachford(y,f3,f1,param);
%
%   These two solvers will converge (up to numerical errors) to the same
%   solution. However, convergence speed might be different. As we perform
%   only 100 iterations with both of them, we do not obtain exactly the
%   same result.
%
%   Figure 5: This figure shows the reconstructed image by solving problem II using the Forward Backward algorithm.
%
%       
%
%   Figure 6: This figure shows the reconstructed image by solving problem II using the Douglas Rachford algorithm.
%
%       
%
%   Remark: The parameter lambda (the regularization parameter) and
%   epsilon (The radius of the l2 ball) can be chosen empirically.
%   Some methods allow to compute those parameter. However, this is far
%   beyond the scope of this tutorial.
%
%   References:
%     P. Combettes and J. Pesquet. Proximal splitting methods in signal
%     processing. Fixed-Point Algorithms for Inverse Problems in Science and
%     Engineering, pages 185-212, 2011.
%     
%
%   Url: http://unlocbox.sourceforge.net/doc/demos/demo_unlocbox.php

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

% Author: Nathanael Perraudin
% Date: fev 23 2012
%

%% Initialisation
clear all;
close all;

% Loading toolbox
init_unlocbox();

verbose=2; % verbosity level

%% Load an image

% Original image
im_original=cameraman; 

% Displaying original image
imagesc_gray(im_original,1,'Original image');  

%% Creation of the problem

sigma_noise = 20/255;
im_noisy=im_original+sigma_noise*randn(size(im_original));

% Create a matrix with randomly 50 % of zeros entry
p=0.5;
matA=rand(size(im_original));
matA=(matA>(1-p));
% Define the operator
A=@(x) matA.*x;

% Depleted image
y=matA.*im_noisy;

% Displaying noiy image
imagesc_gray(im_noisy,2,'Noisy image');  

% Displaying depleted image
imagesc_gray(y,3,'Depleted image');  

%% Setting the proximity operator

% setting the function f1 (norm TV)
paramtv.verbose = verbose-1;
paramtv.maxit = 50;
f1.prox=@(x, T) prox_tv(x, T, paramtv);
f1.eval=@(x) norm_tv(x);   

% setting the function f2 
param_proj.epsilon = sqrt(sigma_noise^2*length(im_original(:))*p);
param_proj.A = A;
param_proj.At = A;
param_proj.y = y;
param_proj.verbose = verbose-1;
f2.prox=@(x,T) proj_b2(x,T,param_proj);
f2.eval=@(x) eps;

% setting the function f3
lambda = 10;
param_l2.A = A;
param_l2.At = A;
param_l2.y = y;
param_l2.verbose = verbose-1;
param_l2.tight = 0;
param_l2.nu = 1;
f3.prox=@(x,T) prox_l2(x,lambda*T,param_l2);
f3.grad=@(x) 2*lambda*A(A(x)-y);
f3.eval=@(x) lambda*norm(A(x)-y,'fro')^2;

%% Solving problem I

% setting different parameters for the simulation
param.verbose = verbose;    % display parameter
param.maxit = 100;    % maximum number of iterations
param.tol = 1e-5;    % tolerance to stop iterating
param.gamma = 1 ;     % Convergence parameter
% solving the problem with Douglas Rachord
sol = douglas_rachford(y,f1,f2,param);

%% Displaying the result
imagesc_gray(sol,4,'Problem I - Douglas Rachford');   

%% Solving problem II (forward backward)
param.gamma=0.5/lambda;     % Convergence parameter
param.tol=1e-5;
% solving the problem with Douglas Rachord
sol21 = forward_backward(y,f1,f3,param);

%% Displaying the result
imagesc_gray(sol21,5,'Problem II - Forward Backward');   

%% Solving problem II (Douglas Rachford)
param.gamma=0.5/lambda;     % Convergence parameter
sol22 = douglas_rachford(y,f3,f1,param);

 %% Displaying the result
imagesc_gray(sol22,6,'Problem II - Douglas Rachford');  

%% Close the UNLcoBoX
close_unlocbox();

