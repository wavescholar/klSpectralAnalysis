%% Color panels by eigenfunction values
pEigenIdxs = 4;
pSamples = 1:5;
psize_Tile = [160,160];

for i = 1:length(pSamples)
    I = imread([DataDir 'MaskData/' Tiles_labels{pSamples(i)} 'Image-MarkerImage.jpg']);
    tile_idxs = find(Tiles_numlabels==pSamples(i));
    fig=figure;
    set(gcf,'Position',[329,132,1178,930]);
    imagesc(I);
    axis off;hold on;
    for k = 1:length(pEigenIdxs),
        %colors=mapvectorstocolors(G.EigenVecs(tile_idxs,pEigenIdxs(k)),colormap);
        colors=mapvectorstocolors(F(tile_idxs),colormap);
        h=scatter(Tiles_loc(1,tile_idxs)*psize_Tile(1),Tiles_loc(2,tile_idxs)*psize_Tile(2),50,colors,'filled');
    end;
    set(gcf,'Name',Tiles_labels{pSamples(i)});
end;