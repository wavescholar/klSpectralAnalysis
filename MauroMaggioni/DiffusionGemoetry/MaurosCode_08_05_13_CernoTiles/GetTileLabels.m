function [TileData_numlabels,TileData_labels] = GetTileLabels( TileData )

for i = 1:length(TileData),
    Labels{i} = TileData{i}{end,1};
end;

[TileData_labels,~,TileData_numlabels] = unique(Labels);

return;