t=[[ 4   460.31  1414.53  2208.14];
   [ 8   211.20   646.59  1064.94];
   [16    97.16   320.78   417.95];
   [32    42.78   137.51   210.84];
   [64    19.50    63.34    83.24];
   [128   13.78    42.71    55.17];
   [256    9.87    31.33    47.95]];

tm = min(min(t(:,2:4)));
tM = max(max(t(:,2:4)));
tm = floor(log2(tm));
tM = ceil(log2(tM));

tl = [tm:tM];
for i = 1:length(tl)
  tl(i) = 2^tl(i);
end

hf = figure;
set(hf,'Position',[500 250 700 850]);
set(hf,'PaperPositionMode','auto');

hp = plot(log2(t(:,1)),log2(t(:,2)),'k');
set(hp,'Marker','.');
hold on
hp = plot(log2(t(:,1)),log2(t(:,3)),'b--');
set(hp,'Marker','.');
hp = plot(log2(t(:,1)),log2(t(:,4)),'r-.');
set(hp,'Marker','.');
axis equal
axis tight
hh = gca;
set(hh,'XLim',[2 8]);
set(hh,'YLim',[tm tM]);
set(hh,'XTickLabel',[4;8;16;32;64;128;256]);
set(hh,'YTickLabel',tl);
xlabel('Number of processors');
ylabel('CPU time (s)');
grid on
box on
legend('STATES','STG','STG\_FULL');




