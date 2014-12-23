function [eigvec eigval] = Diffusion_Family(data,ScaleSigma,knn,flag)

% This function generates the eigen function of the diffusion operator, the
% forward Focker Plank operator and the Laplace Beltrami operator for the
% given data set. The algorithm is the one suggested by coifman and Lafon, 
% 2006.
% Data - N (# of points) x D (dimensions)
% ScaleSigma - Scaling factor for the Gaussian Kernel
% Flag - Df (Diffusion Operator)
%      - Lb (Laplace Beltrami operator, independent of sampling density)
%      - Fp (Forward Focker Plank operator)

% Author - Ambarish Jash

% Creating the adjacency graph with knn number of nearest neighbors.
distances = adjacency(data,'nn',knn);
index     = distances ~= 0;
% Getting average value of distance. Note that those having 0 value are not
% taken into account.
mean_dist = mean(distances (index));
delta     = mean_dist*ScaleSigma;
% Using gaussian kernel
W = zeros(size(distances));
for row_count =1:size(distances,1)
    for col_count = 1:size(distances,2)
        if distances(row_count,col_count) ~=0
            W(row_count,col_count) = exp(-(distances(row_count,col_count)/delta));
        end
    end
end

% Making the weight matrix symmetric
Ws     = (W + W')./2;
if strcmp(flag,'Df')
    for count = 1:size(deg,1)
        W(count,:) = W(count,:)./deg(count,1);
    end
    
    opts.issym  = 1;
    opts.disp   = 0;
    opts.isreal = 1;
    [eigvec eigval] = eigs(W,30,'lm',opts);
    % Getting orthonormal diffusion coordinates
    for count = 1:size(eigvec,2)
        eigvec(:,count) = eigval(count,count)^64.*eigvec(:,count);
        eigvec(:,count) = eigvec(:,count)/norm(eigvec(:,count));
    end
elseif strcmp(flag,'Lb')
    % Eigen Functions of the Laplace Beltrami operator irresoective of
    % sampling density.
    deg    = sum(W,2);
    degree = zeros(size(W));
    for count = 1:size(W,1)
        if deg(count) ~= 0
            degree(count,count) = 1/deg(count);
        end
    end
    clear count
    % In this case alpha is 1
    kernel = degree*Ws*degree;
    deg    = sum(kernel,2);
    for count = 1:size(kernel,1)
        if deg(count) ~= 0
            degree(count,count) = 1/deg(count);
        end
    end
    kernel = degree*kernel*degree;
    
    % Parameteres for eigs
    opts.issym  = 0;
    opts.disp   = 0;
    opts.isreal = 1;
    nEigs       = 50;
    [eigvec, eigval]= eigs (kernel, nEigs,'lm',opts);
    eigvec          = eigvec';
elseif strcmp(flag,'Fp')
    % Eigen function of the Focker Plank operator
    deg    = sum(Ws,2);
    degree = zeros(size(Ws));
    for count = 1:size(Ws,1)
        if deg(count) ~= 0
            degree(count,count) = 1/sqrt(deg(count));
        end
    end
    clear count
    kernel = degree*Ws*degree;
    deg    = sum(kernel,2);
    for count = 1:size(kernel,1)
        if deg(count) ~= 0
            degree(count,count) = 1/sqrt(deg(count));
        end
    end
    kernel = degree*kernel*degree;
    
    %Parameters for eigs
    opts.issym  = 1;
    opts.disp   = 0;
    opts.isreal = 1;
    nEigs       = 50;
    
    [eigvec, eigval] = eigs (kernel, nEigs,'lm',opts);
    eigvec           = eigvec';
    % Normalization from the eigenfunction of the symmeric kernel to the
    % eigen functions of the Focker Plank operator.
    for count = 1:size(eigvec,1)
        eigvec(count,:) = eigvec(count,:)*degree;
        eigvec(count,:) = eigvec(count,:)/norm(eigvec(count,:));
    end;
end
