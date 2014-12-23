% formatscatterplot.m
% -------------------------------------------------------------------------------------------
% formatscatterplot couldn't be easier to use.
% enlarges text and figure outlines
% for current figure.
% run BEFORE adding axis labels, titles,
% etc.
%
% see also: formatheatmap | formatplot | formathist | formatnonnegativeheatmap | formatpdbfig
% by andrej j. savol, dec. 30, 2011 ajs171@pitt.edu
% -------------------------------------------------------------------------------------------

box off;

set(gcf,'Color',[1 1 1]);
p = get(gcf);

% formatplot
for plt = 1 : length(p.Children)
   h = p.Children(plt);
   axes(h)
   % axis square

   set(h,'FontSize',14);
   set(h,'FontName','Helvetica');
   set(h,'FontWeight','Bold');
   set(h,'LineWidth',1);
   set(h,'Color',[1 1 1]);


end

% draw axis lines
for plt = 1 : length(p.Children)
   h = p.Children(plt);

   x = get(h,'XLim');
   y = get(h,'YLim');
   z = get(h,'ZLim');

   X = [x(2) x(2) x(2);
       x(2) x(2) x(1)];

   Y = [y(1) y(1) y(2);
       y(1) y(2) y(2)];

   Z = [z(1) z(2) z(2);
       z(2) z(2) z(2)];

   axes(h)
   line(X,Y,Z,'color','k','LineWidth',1);
end