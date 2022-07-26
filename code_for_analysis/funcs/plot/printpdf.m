function printpdf(h,outfilename,res)
% 
% printpdf(gcf,'trash','-r600')
% 

set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos = get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'auto');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename,res);