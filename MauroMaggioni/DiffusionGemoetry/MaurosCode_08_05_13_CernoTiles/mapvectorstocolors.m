function C = mapvectorstocolors( V,colors )

uV = unique(V);
if length(uV)<=2,
    idxs(V<=0)=1;
    idxs(V>0)=size(colors,1);
    C = colors(idxs,:);
else
    C = colors(round((V-min(V))/(max(V)-min(V))*(size(colors,1)-1))+1,:);
end;

return;