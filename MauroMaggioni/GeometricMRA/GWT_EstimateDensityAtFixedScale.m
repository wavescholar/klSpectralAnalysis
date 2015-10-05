function DensEst = GWT_EstimateDensityAtFixedScale( GWT, Data, j, opts )

if nargin<4,
    opts =[];
end;

if ~isfield(opts,'MaxPtsForDensityEst'),    opts.MaxPtsForDensityEst = inf; end;

%% Find the planes at scale j
cp_idx = get_partition_at_scale(GWT,j);
n_nodes = length(cp_idx);

%% Construct a density approximation at scale j, based on TrainingData on each plane
% Allocate memory
DensEst.Density(n_nodes)    = kde;
DensEst.Planes              = cell(1,n_nodes);
%ent                         = zeros(1,n_nodes);

for k = 1:n_nodes
    % Density estimation for the scaling function coefficients
    % Find the scaling coefficients points belonging to each node at scale j
    [lScalCoeffs] = GWT_getScalCoeffs_atnode  ( GWT, Data, cp_idx(k) );
    % Estimate empirical probability of how many points belong to each plane
    DensEst.pi(k)   = size(lScalCoeffs,2);
    % Model the distribution of scaling coefficients
    if DensEst.pi(k)>0,
        if opts.MaxPtsForDensityEst == inf,
            DensEst.Density(k)          = kde(lScalCoeffs,'rot');                                   % MM_DBG: figure;plot(lScalCoeffs(1,:),lScalCoeffs(2,:),'.');samples=sample(DensEst,1000);hold on;plot(samples(1,:),samples(2,:),'r.');
        else
            tmp                         = randperm(size(lScalCoeffs,2));
            DensEst.Density(k)          = kde(lScalCoeffs(:,tmp(1:min([size(lScalCoeffs,2),opts.MaxPtsForDensityEst]))),'rot');
        end;
        DensEst.Planes{k}           = GWT.ScalFuns{cp_idx(k)};
        DensEst.PlaneCenter(:,k)    = GWT.Centers{cp_idx(k)};
        DensEst.cp_idx(k)           = cp_idx(k);
        %ent(k) = entropy(DensEst.Density(k));
    end;
    % Find the scaling coefficients points belonging to each node at scale j
    [lWavCoeffs,lWav_cp_idx] = GWT_getWavCoeffs_atnode   ( GWT, Data, cp_idx(k) );
    % Estimate empirical probability of how many points belong to each plane
    for r = 1:length(lWavCoeffs),
        DensEst.pi_w(k,r)   = size(lWavCoeffs{r},2);
        % Model the distribution of scaling coefficients
        if DensEst.pi_w(k,r)>0,
            if opts.MaxPtsForDensityEst == inf,
                DensEst.Density_w(k,r)      = kde(lWavCoeffs{r},'rot');                              % MM_DBG: figure;plot(lScalCoeffs(1,:),lScalCoeffs(2,:),'.');samples=sample(DensEst,1000);hold on;plot(samples(1,:),samples(2,:),'r.');
            else
                tmp                         = randperm(size(lWavCoeffs{r},2));
                DensEst.Density_w(k,r)      = kde(lWavCoeffs{r}(:,tmp(1:min([size(lWavCoeffs{r},2),opts.MaxPtsForDensityEst]))),'rot');
            end;
            DensEst.Planes_w{k,r}           = GWT.WavBases{lWav_cp_idx(r)};
            DensEst.PlaneCenter_w(:,k,r)    = GWT.Centers{lWav_cp_idx(r)}+GWT.WavConsts{lWav_cp_idx(r)};
            DensEst.cp_idx_w(k,r)           = lWav_cp_idx(r);
        end;
    end;
end;

%% Renormalize pi
DensEst.pi = DensEst.pi/sum(DensEst.pi);

%% Compute (differential) entropy of distribution
%DensEst.entropy = entropy(DensEst.pi) + sum(DensEst.pi.*ent);

%% Return
return;