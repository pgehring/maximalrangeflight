clear;
clc;
close all;

%IP = readmatrix('test_4_2_IP.txt');
SQP = readmatrix('./results/test_0_1_SQP.txt');
IP = readmatrix('./results/test_0_1_IP.txt');

set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLineMarkerSize',4);
set(0, 'DefaultLineLineWidth', 1);
FigW=16;
FigH=8;
f = figure(1);
set(f,'defaulttextinterpreter','latex',...
            'PaperUnits','centimeters','PaperSize',[FigW FigH],...
                'PaperPosition',[0,0,FigW,FigH],'Units','centimeters',...
                'Position',[0,0,FigW,FigH]);

set(f, 'PaperPositionMode', 'manual');

subplot(2,1,1);
plot(SQP(:,1),SQP(:,3),"-r");
hold on
plot(IP(:,1),IP(:,3),"-b");
%title('Zielfunktion $F$')
ylabel('Zielfunktion $F$ in $[m]$')
xlabel("Anzahl der Iterationen")


subplot(2,1,2);
semilogy(SQP(:,1),SQP(:,4),"-r");
hold on 
semilogy(IP(:,1),IP(:,4),"-b");
%title('Feasibility')
ylabel('Zul{\"a}ssigkeit in [1]')
xlabel("Anzahl der Iterationen")


% legend1 = legend(gca,["SQP", "IP"]);
% % set the legend where I really want it.
% set(legend1,'Units','centimeters','Position',[8 3 1 1],'Orientation','horizontal');
      


ax = findobj(f, 'type', 'axes');

%ax(1).YAxis.Exponent = 0;

set(f, 'PaperPositionMode', 'auto');

ax(1).Units = 'centimeters';
test=ax(1).Position;
ax(1).Position=[test(1),test(2),test(3),test(4)];
set(findall(gcf,'-property','FontSize'),'FontSize',8);

legend(ax(1),["SQP", "IP"], 'Interpreter','latex','Location','southeast','FontSize',6,'Orientation','horizontal');
legend(ax(2),["SQP", "IP"], 'Interpreter','latex','Location','southeast','FontSize',6,'Orientation','horizontal');
print(f,'-dpdf','-r600','./results/Vergleich_SQP_IP.pdf');
