% formatfig.m
% --------------------------------------
% enlarges text and figure outlines
% for current figure. 
% run BEFORE adding axis labels, titles,
% etc.
% 
% see also: formatheatmap
% --------------------------------------

set(gca,'FontSize',24);
set(gca,'FontName','Helvetica');
set(gca,'FontWeight','Bold');
set(gca,'LineWidth',2);
set(gcf,'Color',[1 1 1]);
box on