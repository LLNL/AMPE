%%% Plot parallel timing results for brusselator example
%%% (simulatino and ASA)
%%%

statsSIM = load('bruss.sim');
colSIM = 4;   %% column with CPU times

statsASA = load('bruss.asa');
colASA = size(statsASA,2);  %% column with CPU times


procsSIM = unique(statsSIM(:,1));
for i = 1:length(procsSIM)
  np = procsSIM(i);
  nt = statsSIM(find(statsSIM(:,1)==np),colSIM);
  cpuSIM(i) = mean(nt);
end

procsASA = unique(statsASA(:,1));
for i = 1:length(procsASA)
  np = procsASA(i);
  nt = statsASA(find(statsASA(:,1)==np),colASA);
  cpuASA(i) = mean(nt);
end

cpuSIM_max = max(statsSIM(:,colSIM));
cpuASA_max = max(statsASA(:,colASA));
cpu_max = max(cpuSIM_max,cpuASA_max);

figure
set(gcf,'position',[500 400 950 500]);
hold on


plot(statsSIM(:,1),statsSIM(:,colSIM),'bd');
plot(statsASA(:,1),statsASA(:,colASA),'rd');


procs = unique([procsSIM ; procsASA]);

set(gca,'XTick',procs);
set(gca,'XLim',[procs(1)-1 procs(end)+1]);
set(gca,'YLim',[0 1.2*cpu_max]);

set(gca,'FontName','times')
%%set(gca,'FontWeight','bold')
set(gca,'FontSize',20)

set(gca,'Color',[0.85 0.85 1.0 ]);

%% Some nicer grid lines

xtick = get(gca,'XTick');
ytick = get(gca,'YTick');

xlim = get(gca,'XLim');
ylim = get(gca,'YLim');

for i = 1:length(xtick)
  plot([xtick(i) xtick(i)],ylim,'color',[0.7 0.7 0.7]);
end
for i = 1:length(ytick)
  plot(xlim, [ytick(i) ytick(i)],'color',[0.7 0.7 0.7]);
end

%% The box
plot([xlim(1) xlim(1)],ylim,'k');
plot([xlim(2) xlim(2)],ylim,'k');
plot(xlim,[ylim(1) ylim(1)],'k');
plot(xlim,[ylim(2) ylim(2)],'k');



hp1 = plot(procsSIM, cpuSIM,'bd-','LineWidth',2);
hp2 = plot(procsASA, cpuASA,'rd-','LineWidth',2);


hl = legend([hp1,hp2],'Simulation','ASA');
set(hl,'color','white');
xlabel('Number of processes');
ylabel('CPU time (s)');

