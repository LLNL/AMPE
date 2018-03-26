load cvsdiscx.mat
figure
hp=plot(s1(:,1),s1(:,2));
set(hp,'LineWidth',2)
hx=xlabel('time');
box on
set(gca,'FontSize',14);
set(hx,'FontSize',14);

figure
hold on
hp1=plot(s2(:,1),s2(:,2));
set(hp1,'LineWidth',2)
hp2=plot(s3(:,1),s3(:,2),'r');
set(hp2,'LineWidth',2)
box on
hx=xlabel('time');
set(gca,'FontSize',14);
set(hx,'FontSize',14);
   
magnify

