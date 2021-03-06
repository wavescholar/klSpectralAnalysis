function [lambda_out,uu2] = eigenfunctions_image(DATA,SIGMA,NUM_EVECS,EPSILON,figs)
  %
  %   [lambda_out,uu2] = eigenfunctions(DATA,SIGMA,NUM_EVECS,EPSILON)  
  %  
  % Function that computes approximate smallest NUM_EVECS of graph
  % Laplacian of DATA, using eigenfunctions approach from NIPS09 paper:
  % "Semi-supervised learning in gigantic image collections" 
  % by R. Fergus, Y. Weiss and A. Torralba. Note that this method assumes
  % that the distribution of DATA is seperable.
  %
  % Inputs: 
  %   1. DATA - nPoints x nDims real matrix. 
  %   2. SIGMA - 1 x 1 real. Affinity between points in graph Laplacian
  %   3. NUM_EVECS -  1 x 1 integer. Number of approximate eigenvectors
  %   to compute.
  %   4. EPSILON (optional) - 1 x 1 real. Default = 0.1. Offset added to
  %   histogram in numericalEigenFunctions to avoid divide by 0 problems.
  %
  %
  % Outputs:
  %   1. lambda_out - NUM_EVECS x NUM_EVECS diagonal matrix holding
  %   approximate eigenvalues of smallest NUM_EVECS eigenvectors.
  %   2. uu2 - nPoints x NUM_EVECS real matrix of approximate
  %   eigenvectors. 
  %
  %  Example usage:
  %   % Compute smallest 5 approximate eigenvectors of graph Laplacian,
  %   % using affinity scaling of 0.2:
  %   [d,u] = eigenfunctions(data,0.2,5);  
  %  
  % This implementation of the numerical eigenfunctions approach is general
  % and can be applied to real data (it is a cleaned up version of the
  % code used to generate the results in the NIPS paper). Note that it
  % does assume seperability of the data, so you might want to preprocess
  % it first with PCA. Note that whatever transformation you use must not
  % change the distances between points, thus only rotations are
  % permitted and things like ICA are not allowed.
  %
  % If you have new data for which you want to use the same set of
  % eigenfunctions, then you just need to do the clipping and the
  % interpolation in the existing eigenfunctions (stages 3 & 4) below.  
  %  
  % If you do use this code in your paper, please cite our paper. For your
  % convenience, here is the bibtex: 
  %
  % @inproceedings{FergusWeissTorralbaNIPS09,
  %  author = "Fergus, R. and Weiss, Y. and Torralba, A.",
  %  title = "Semi-supervised Learning in Gigantic Image Collections",
  %  booktitle = "NIPS",
  %  year = "2009"
  % } 
  %
  % Version 1.0. Rob Fergus & Yair Weiss. 11/30/09.  
   
    
%   % To avoid crazy values for points on extreme end of each 1D distribution we
%   % clip eigenvalues of points in the bottom and top CLIP_MARGIN percent
%   % of points. This is important for stability of the semi-supervised
%   % learning solution. 
%   CLIP_MARGIN = 2.5; %%% percent  
  
if nargin < 5
    figs = 0; % don't plot figures!
end

% working with an image:
  nPoints = prod(size(DATA));
  nDims = ndims(DATA);
    
  if nargin==3
    %%% default value for numericalEigenFunctions
    EPSILON = 0.1;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% 1. Compute eigenfunctions seperately for each input dimension
  
  fprintf('Computing eigenfunctions\n');
  %%% compute numerical eigenfunctions per dimension
  Cbins = cell(nDims,1);
  Clambdas = cell(nDims,1);
  Cpp = cell(nDims,1);
  for a=1:nDims
    % Solve eqn. 2 in the NIPS paper, one dimension at a time.
    dimmm = sum(DATA,a);
    figure; plot(dimmm); title(['sums across dimension ',num2str(a)],'FontSize',20);
    [Cbins{a},Cg{a},Clambdas{a},Cpp{a}] = numericalEigenFunctions_image(dimmm,SIGMA,EPSILON,figs);
    %[bins(:,a),g(:,:,a),lambdas(:,a),pp]=numericalEigenFunctions(DATA(:,a),SIGMA,EPSILON);
  end
%   while size(Cbins{1},1) ~= size(Cbins{2},1)
%       nDims = min(size(Cbins{1},1),size(Cbins{2},1));
%       for a=1:nDims
%           % Solve eqn. 2 in the NIPS paper, one dimension at a time.
%           [Cbins{a},Cg{a},Clambdas{a},Cpp{a}] = numericalEigenFunctions_vmb(DATA(:,a),SIGMA,EPSILON,ndims);
%           %[bins(:,a),g(:,:,a),lambdas(:,a),pp]=numericalEigenFunctions(DATA(:,a),SIGMA,EPSILON);
%       end
%       
%   end
  nbins = size(Cbins{1},1);
  bins = zeros(nbins,nDims);
  g = zeros(nbins,nbins,nDims);
  lambdas = zeros(nbins,nDims);
  for a = 1:nDims
      bins(:,a) = Cbins{a};
      g(:,:,a) = Cg{a};
      lambdas(:,a) = Clambdas{a};
  end
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% 2. Choose smallest NUM_EVECS eigenfunctions over all input dimensions
  
  %%% now get all eignvalues and find smallest k
  all_l = lambdas(:);
  
% vmb: 
%   %%% eliminate all the really tiny ones by setting them to be huge
%   q=find(all_l<1e-10);
%   all_l(q) = 1e10;
  
  %%% sort for smallest 
  [lambda_out,ind] = sort(all_l);
  
  % vmb:
  % if both directions have a zero eigenvalue, only keep smallest:
  mindims = min(lambdas)
  if size(find(mindims < 10^-12),2) == nDims % if both minima are practically zero
      disp('minimum eigval in both dimensions is 0, so only keeping smaller 0');
     % remove second value from lambda_out
     lambda_out = lambda_out([1 3:end]);
     ind = ind([1 3:end]);
  end
     
  
  
  %%% take first NUM_EVECS
  use_l = ind(1:NUM_EVECS);
  
  %%% get indices of picked eigenvectors
  lambda_out = diag(lambda_out(1:NUM_EVECS));
  [ii,jj]=ind2sub(size(lambdas),use_l);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% 3. Do clipping of data distribution per dimension
%   LOWER = CLIP_MARGIN/100; %% clip lowest CLIP_MARIN percent
%   UPPER = 1-LOWER; % clip symmetrically
%   
%   for a=1:nDims
%     [clip_lower(a),clip_upper(a)] = percentile(DATA(:,a),LOWER,UPPER);
%     q = find(DATA(:,a)<clip_lower(a));
%     % set all values below threshold to be constant
%     DATA(q,a) = clip_lower(a);
%     q2 = find(DATA(:,a)>clip_upper(a));
%     % set all values above threshold to be constant
%     DATA(q2,a) = clip_upper(a);
%   end

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% 4. Interpolate DATA in smallest NUM_EVECS eigenfunctions
  
  %%% define memory first, single to save memory
  uu2 = zeros(nPoints,NUM_EVECS,'single');
    
  fprintf('Interpolating data in eigenfunctions\n');
  %%%% Now do the interpolation
  for a=1:NUM_EVECS

    %% get bin indices
    bins_out(:,a) = bins(:,jj(a));
    %% and function value
    uu(:,a) = g(:,ii(a),jj(a)); % ii is for dim 1, jj is for dim 2
 
    %%% now look up data in it.
    % get all points from this bin - bin(ii(a)) from dim(jj(a))
    if jj(a) == 1
        dimmm = DATA'; % rastorize row-wise
    elseif jj(a) == 2 
        dimmm = DATA; % rastorize column-wise
    end
 %   uu2(:,a) = interp1(bins_out(:,a),uu(:,a),DATA(:,jj(a)),'linear','extrap');
 %   uu2(:,a) = interp1(bins_out(:,a),uu(:,a),dimmm(:),'linear','extrap');
    uu2(:,a) = interp1(bins_out(:,a),uu(:,a),dimmm(:),'linear','extrap');
     %%% matlab implementation is rather slow. 
  
  end
    

  %%% make sure approx. eigenvectors are of unit length.
  suu2 = sqrt(sum(uu2.^2));
  uu2 = uu2 ./ (ones(nPoints,1) * suu2);

