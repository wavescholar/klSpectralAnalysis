function h2 = figurev(id)

if nargin < 1
    h2 = figure;
else
    h2 = figure(id);
end
% set(h2,'Position',[ 127   753   861   717]); % left screen
%set(h2,'Position',[ 2000   600   861   717]); % right screen
%set(h2,'Position',[ 1701 146  1615  804]); % right screen
%set(h2,'Position',[ 48 106  1186  574]); % macbook air
set(h2,'Position',[14 47 1779 912]); % asus
%set(gca,'FontSize',12);

% set(gca,'FontSize',24);
% set(gca,'FontName','Helvetica');
% set(gca,'FontWeight','Bold');
% set(gca,'LineWidth',2);
% set(gcf,'Color',[1 1 1]);
% box on