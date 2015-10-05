function [risk,risk_labels,diagnosis,diagnosis_labels] = GetRiskAndDiagnosis( TileData_labels )

for i = length(TileData_labels):-1:1,
    idxs = find(TileData_labels{i}=='-');
    risk_alllabels{i}       = TileData_labels{i}(1:idxs(1)-1);
    diagnosis_alllabels{i}  = TileData_labels{i}(idxs(1)+1:idxs(2)-1);
end;

[risk_labels,~,risk]            = unique(risk_alllabels);
[diagnosis_labels,~,diagnosis]  = unique(diagnosis_alllabels);


return;