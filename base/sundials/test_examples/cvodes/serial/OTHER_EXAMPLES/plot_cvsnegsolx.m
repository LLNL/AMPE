load('no_cntrl.out');
load('cntrl.out');

figure;
hold on
plot([1e9 4e10],[0 0],'k');
hp1=plot(no_cntrl(:,1),no_cntrl(:,2),'r');
set(hp1,'LineWidth',2);
hp2=plot(cntrl(:,1),cntrl(:,2),'b');
set(hp2,'LineWidth',2);

set(gca,'XScale','log');
set(gca,'XLim',[1e9,4e10]);
set(gca,'YLim',[-0.5e-6,1e-6]);
set(gca,'FontSize',14);
grid on
box on
hl = xlabel('time');
set(hl,'FontSize',14);
hl = ylabel('y_1');
set(hl,'FontSize',14);

legend([hp1,hp2],'no control','with control');