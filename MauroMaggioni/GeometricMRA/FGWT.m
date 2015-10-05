function Data = FGWT(gW, x)

%
% Fast Geometric Wavelets Transform
%
% Input:
%     gW: the geometric wavelets structure computed by the function
%         geometric_wavelets_transformatin. In fact, only the following fields
%         are needed:
%       .cp: the vector encoding the metis tree structure
%       .Scales: vector of scales of the nets
%       .Centers: the local centers of the nets
%       .WavBases: wavelet bases
%       .WavConsts: translations associated with the wavelet bases
%   x: a new point, 1-by-D vector
% 
% Output: 
%   Data: structure of the following fields:
%       .WavCoeffs: cell array of wavelet coefficients
%       .chain: path from the root to the leaf node along the tree 
%               This auxillary field can be used to extract from gW 
%               the sequences of wavelet bases and translations by
%               using gW.WavBases{chain} and gW.WavConsts{chain}

J = max(gW.Scales); % number of scales

Data            = struct();
Data.ScalCoeffs = cell(1,J); 
Data.WavCoeffs  = cell(1,J); 

%% Find the leaf node that is closest to x
net = find_nearest_leaf_node(gW, x); 

%% Compute transform bottom up
j = gW.Scales(net);

Data.chain = zeros(1,j);

if  j==1 % only one scale
    
    Data.chain(1)       = net;
    Data.WavCoeffs{1}   = (x-gW.Centers{net}) * gW.WavBases{net};
    Data.ScalCoeffs{1}  = Data.WavCoeffs{1};
    
else    % go bnottom up from the leaf to the root of the GMRA tree
    
    iFineNet    = net; % current scale
    iCoarseNet  = gW.cp(net); % next scale
    
    if ~gW.opts.orthogonalizing
        if gW.opts.avoidLeafnodePhi
            finestBasis = [gW.ScalFuns{iCoarseNet} gW.WavBases{iFineNet}];
        else
            finestBasis = gW.ScalFuns{iFineNet};
        end
        % x_J: finest approximation at scale J
        % x_j: scale j approximation for any j
        Data.ScalCoeffs{j} = (x-gW.Centers{iFineNet})*finestBasis;
        x_J = Data.ScalCoeffs{j}*finestBasis' + gW.Centers{iFineNet}; % best approximation
        x_j = x_J;
    else
        Data.ScalCoeffs{j} = (x-gW.Centers{iFineNet})*gW.ScalFuns{iFineNet};
    end
    
    while j>1

        Data.chain(j) = iFineNet;
        
        if ~gW.opts.orthogonalizing
            wavelet = gW.WavConsts{iFineNet};
            
            if ~isempty(gW.WavBases{iFineNet})
                Data.WavCoeffs{j} = (x_j-gW.Centers{iFineNet}) * gW.WavBases{iFineNet};
                wavelet = wavelet + Data.WavCoeffs{j}*gW.WavBases{iFineNet}';
            end
            
            if gW.opts.addTangentialCorrections
                wavelet = wavelet - (x_J-x_j)*gW.ScalFuns{iCoarseNet}*gW.ScalFuns{iCoarseNet}';
            end
            
            x_j = x_j - wavelet;
            Data.ScalCoeffs{j-1} = (x_j-gW.Centers{iCoarseNet})*gW.ScalFuns{iCoarseNet};
        else
            Data.ScalCoeffs{j-1} = (x-gW.Centers{iCoarseNet})*gW.ScalFuns{iCoarseNet};
            Data.WavCoeffs{j} = (x-gW.Centers{iFineNet}) * gW.WavBases{iFineNet};   
        end
        
        j = j-1;
        iFineNet = iCoarseNet;
        iCoarseNet = gW.cp(iFineNet);
        
    end
    
    % j = 1
    Data.chain(1) = iFineNet;
    Data.WavCoeffs{1} = (x_j-gW.Centers{iFineNet}) * gW.WavBases{iFineNet};
    
end

return;
