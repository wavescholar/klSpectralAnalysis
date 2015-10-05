function paintArvind(snapEnergy,data, ploton);
% paints colors in snapEnergy onto data

if nargin < 3
    ploton = 0;
end

[min(snapEnergy) max(snapEnergy)]

%[sssE, ssiE] = sort(snapEnergy);
noOfBins = 21;
[nn,ss] = hist(snapEnergy,noOfBins);

if ploton
    figure; subplot(2,1,1); plot(ss,nn/sum(nn));
    subplot(2,1,2); semilogy(ss,nn/sum(nn)); grid on
end

% thresholding the tails: include only those that are within 25% of the
% center 
diff = max(snapEnergy) - min(snapEnergy);
thr = .25*diff;
snapEnergy(snapEnergy < min(snapEnergy) + thr) = min(snapEnergy) + thr;
snapEnergy(snapEnergy > max(snapEnergy) - thr) = max(snapEnergy) - thr;
diff = max(snapEnergy) - min(snapEnergy);
if diff < eps
    noOfBins = 2;
else
    noOfBins = 101;
end

[nn,ss] = hist(snapEnergy,noOfBins);
% bin width
binWidth = ss(2)-ss(1);
% from bin center to bin edge
ss = ss - binWidth;
% find out which bin is snapEnergy(i) falling in..
snapEnergyBin = floor((snapEnergy - ss(1))/binWidth) + 1;
snapEnergyBin(snapEnergyBin>noOfBins) = noOfBins;


%% plot:

cmap = jet(noOfBins);

pp = [1:length(data)];%[1:20200];%structIds{iter,1};
% snapEnergy = sum(totalNew(:,pp), 1); % not needed since the totalEnergy
% is itself all energy data.

%figure;
qq = data;
for i = 1:length(qq)
    %idx = ssiE(i);
    %h = plot3(qq(ki{iter}(1), idx), qq(ki{iter}(2), idx), qq(ki{iter}(3), idx), 'o', 'color', cmap(i, :), 'MarkerSize', 6, 'MarkerFaceColor', cmap(i, :));
    
    if size(data,2) > 2
        plot3(qq(i,1), qq(i,2), qq(i,3), '.', 'color', cmap(snapEnergyBin(pp(i)), :));%, 'MarkerFaceColor', cmap(snapEnergyBin(pp(i)), :));
    else
        plot(qq(i,1), qq(i,2), '.', 'color', cmap(snapEnergyBin(pp(i)), :));% , 'MarkerFaceColor', cmap(snapEnergyBin(pp(i)), :));
    end
    hold on;
end

% set(gca, 'FontSize', 24);
% set(gca, 'FontWeight', 'bold');
% set(gca, 'LineWidth', 2.0);
% set(gca, 'Box', 'on');
%title(sprintf('level: %d',iter));
% grid on;
axis equal;
axis off; 
%axis tight;
%   saveas(h, sprintf('data/2KKJ_01/ACEMD-prep/ANAL/dPCA_2KKJ_radgyr_%d.fig', iter), 'fig');
%   saveas(h, sprintf('data/2KKJ_01/ACEMD-prep/ANAL/dPCA_2KKJ_radygyr_%d.ai', iter), 'ai');
%   saveas(h, sprintf('data/2KKJ_01/ACEMD-prep/ANAL/dPCA_2KKJ_radgyr_%d.png', iter), 'png');
