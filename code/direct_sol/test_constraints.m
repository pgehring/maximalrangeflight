% test_constraints.m
% Description:
%   Test if a state reaches the maximum back pressure condition of the
%   direct solution.
% Date:         27.08.2021
% Author:       Gehring, Philipp / Karus, Heiko / Goetz, Felix

clear  variables; 
close  all;
clc;

%% Load data
test_results_name = 'test_1_1';
run(strcat('./config/',test_results_name));
data = readmatrix(strcat('./results/',test_results_name,'.txt'));

%% Test if a state reaches the maximum back pressure condition
qmax_q_diff = zeros(N,3);
nof_exc = 0;
for i = 1:N
    q = (prob.alpha*exp(-prob.beta*data(i,1))*data(i,4)^2)/2;
    qmax_q_diff(i,:) = [prob.q_max,q,prob.q_max-q];
    if (q > prob.q_max)
        nof_exc = nof_exc + 1;
    end
end

%% Plot
fig = figure(1);
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');


set(0,'DefaultLineMarkerSize',4);
set(0, 'DefaultLineLineWidth', 1);
FigW=16;
FigH=7;





x = linspace(0,N-1,N);
plot(x,qmax_q_diff(:,1),'r',x,qmax_q_diff(:,2),'b--');
xlabel('Diskretisierungspunkte');
ylabel('$q(v(t),h(t))$');
set(findall(gcf,'-property','FontSize'),'FontSize',8);

legend('$q_{\max}$','$q(v(t),h(t))$','Interpreter','latex','Location','southoutside','FontSize',7,'Orientation','horizontal');
title(['Staudruck \"Ueberpr\"ufung (',num2str(nof_exc),' \"Uberschreitungen)'],'Interpreter','latex','FontSize',10);

set(fig,'defaulttextinterpreter','latex',...
            'PaperUnits','centimeters','PaperSize',[FigW FigH],...
                'PaperPosition',[0,0,FigW,FigH],'Units','centimeters',...
                'Position',[0,0,FigW,FigH]);
set(fig, 'PaperPositionMode', 'auto');

 ax = findobj(fig, 'type', 'axes');
    ax.Units = 'centimeters';
    test=ax.Position;
    ax.Position=[test(1),test(2)-0.35,test(3),test(4)];

print(fig,'-dpdf','-r600',strcat('./results/',test_results_name,'_staudruck.pdf'));
% savefig(fig,strcat('./results/',test_results_name,'_staudruck.fig'));
% saveas(fig,strcat('./results/',test_results_name,'_staudruck.png'));