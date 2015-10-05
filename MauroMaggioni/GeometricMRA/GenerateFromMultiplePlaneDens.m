function [X] = GenerateFromMultiplePlaneDens( MultiplePlanesDens, n, opts )

if nargin<3,                        opts = [];         end;
if ~isfield(opts,'alpha'),          opts.alpha =  1;   end;

X = zeros(size(MultiplePlanesDens.Planes{1},1),n);

n_planes = length(MultiplePlanesDens.pi);

% To speed this up, we first decide how many points are going to be sampled from each plane
plane_assgn = sample(MultiplePlanesDens.pi,n);%mnrnd(1,MultiplePlanesDens.pi,n)+1;
% Find which points are in each plane
for k = n_planes:-1:1,
    lPlaneIdxs{k} = find(plane_assgn==k);
end;

% Now draw points from each plane
for k = 1:n_planes,
    if ~isempty(lPlaneIdxs{k}),
        % Draw from the distribution of the k-th plane
        X(:,lPlaneIdxs{k}) = bsxfun(@plus,opts.alpha*MultiplePlanesDens.Planes{k}*sample(MultiplePlanesDens.Density(k),length(lPlaneIdxs{k})),MultiplePlanesDens.PlaneCenter(:,k));
%         clf;plot3(MultiplePlanesDens.PlaneCenter(1,:),MultiplePlanesDens.PlaneCenter(2,:),MultiplePlanesDens.PlaneCenter(3,:),'.');hold on;plot3(MultiplePlanesDens.PlaneCenter(1,k),MultiplePlanesDens.PlaneCenter(2,k),MultiplePlanesDens.PlaneCenter(3,k),'ro');plot3(X(1,lPlaneIdxs{k}),X(2,lPlaneIdxs{k}),X(3,lPlaneIdxs{k}),'g+');
%         set(gca,'CameraPosition',[0.75505 -8.14674 -0.444029]);title(sprintf('%d',k));pause;
%         for l=1:length(lPlaneIdxs{k}),
%             if norm(X(:,lPlaneIdxs{k}(l))-[0;0;0.5;zeros(7,1)])<1e-1,
%                 k,
%             end;
%         end;
    end;
end;

for l = 1:size(X,2),
    if norm(X(:,l))==0,
        break;
    end;
end;

return;