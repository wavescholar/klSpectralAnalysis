function TileLoc = GetTileLocations( TileData )

TileLoc = zeros(2,length(TileData));

progressbar;
for i = 1:length(TileData),
    TileLoc(:,i) = TileData{i}{118,1}';
    progressbar(i/length(TileData));
end;

return;