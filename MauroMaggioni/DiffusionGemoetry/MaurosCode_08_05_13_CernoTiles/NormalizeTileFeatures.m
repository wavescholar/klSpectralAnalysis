function TileData_norm = NormalizeTileFeatures( TileData, Opts )

%
% function TileData_norm = NormalizeTileFeatures( TileData )
%
% IN:
%   TileData    : TileData for a panel, in Cernostics format
%   [Opts]      : structure of options.
%                   [NormalizationType] : any of {'mean'}. Default: 'mean'.
%
% OUT:
%   TileData_norm : a matrix of normalized TileData
%
%

if nargin<2,    Opts = []; end;
if ~isfield(Opts,'NormalizationType')   Opts.NormalizationType = 'mean'; end;

switch( lower(Opts.NormalizationType) )
    case 'mean'
        % Do the simplest thing possible: replace each vector in the tile data by its mean
        TileData_norm = zeros(length(TileData),117);
        %progressbar;
        parfor i = 1:length(TileData)
            for k = 1:117
                TileData_norm(i,k) = mean(TileData{i}{k,1});
            end
           % progressbar(i/length(TileData));
        end
        
        % Now put all the quantities on the same scale, by zscore
        TileData_norm = zscore(TileData_norm);
        
        TileData_norm = TileData_norm';
    otherwise
        error(sprintf('NormalizationType %s unknown',Opts.NormalizationType));
end;


return;