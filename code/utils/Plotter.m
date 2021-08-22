% Plotter.m:

% Klasse zur graphischen Dartstellung der Ergebnisse

% Date:         27.08.2021
% Author:       Gehring, Philipp / Karus, Heiko / Götz, Felix

classdef Plotter
    properties
        nfigures = 0;
        figures;
        axes;
        n;
        m;
    end
    methods
        %% Konstruktor
        function obj = Plotter()
            close("all");
            set(0,'defaulttextinterpreter','latex');
            set(0,'defaultAxesTickLabelInterpreter','latex');
            set(0,'DefaultLineMarkerSize',4);
            set(0, 'DefaultLineLineWidth', 1);
        end
        
        %%
        function fig = plot_fmincon(obj,t,sol,results_name,titles,labels,order,frame_prop,line_style)
            obj.nfigures = obj.nfigures +1;
            fig = create_fig(obj, obj.nfigures);
            
            FigW=16;
            FigH=11;
            set(fig,'defaulttextinterpreter','latex','PaperUnits','centimeters','PaperSize',[FigW FigH],...
                'PaperPosition',[0,0,FigW,FigH],'Units','centimeters',...
                'Position',[0,0,FigW,FigH]);
            set(fig, 'PaperPositionMode', 'auto');
            
            fig_title_str = strcat('Versuch: \, ','\verb|',results_name,'|');
            fig_title = sgtitle(fig,fig_title_str,'FontSize',12);
            fig_title.Interpreter = 'latex';
            
            [r,c] = size(sol);
            for i = [1:c]
                sol_id = order(i);
                sp = create_subplot(obj,t,sol(:,sol_id),[2, 3],i,titles(sol_id),labels(sol_id),frame_prop(sol_id),line_style(sol_id));
            end
            
            subplots = findobj(fig,'type','axes');
            for i = [1:c]
                sp = subplots(i);
                sp.Units = 'centimeters';
                test=sp.Position;
                if i == 1
                    sp.Position=[test(1)+0.7,test(2)-0.3,test(3),test(4)]; % Steuerung 2
                    ylim(sp,[0-1.48*0.05,1.48+1.48*0.05]);
                    sp.YTick=0:0.37:1.48;
                elseif i == 2
                    sp.Position=[test(1),test(2)-0.3,test(3),test(4)]; % Geschwindigkeit
                elseif i == 3
                    sp.Position=[test(1)-0.7,test(2)-0.3,test(3),test(4)]; % Anstellwinkel
                elseif i == 4
                    sp.Position=[test(1)+0.7,test(2)-0.5,test(3),test(4)]; % Steuerung 1
                    ylim(sp,[0-1260000*0.05,1260000+1260000*0.05]);
                    sp.YTick=0:180000:1260000;
                elseif i == 5
                    sp.Position=[test(1),test(2)-0.5,test(3),test(4)]; % Flughöhe
                else
                    sp.Position=[test(1)-0.7,test(2)-0.5,test(3),test(4)]; % Strecke
                end
            end
        end
        
        function fig = plot_state(obj, t, X, titles, labels)
            % function to plot the state vector over time
            obj.nfigures = obj.nfigures +1;
            fig = create_fig(obj, obj.nfigures);
            %
            [r,c] = size(X);
            for i = [1:c]
                sp = create_subplot(obj,t,X(:,i),[c, 1],i,titles(i),labels(i),0.5,"b-");
            end
        end
        
        function fig = create_fig(obj, fig_id)
            fig = figure();
            obj.figures(fig_id) = fig;            
        end
        
        function ax = create_subplot(obj,t,y,pos,ax_id,fig_title,fig_label,frame_prop,line_style)
            ax = subplot(pos(1), pos(2), ax_id);
            obj.axes(ax_id) = ax;
            if (nargin <= 6)
                plot(t, y, '-k')
            else
                plot(t, y, line_style)
            end
            y_min = min(y);
            y_max = max(y);
            p = 0.05;
            if abs(y_min) < abs(y_max)
                y_min = y_min - y_max*p;
                y_max = y_max + y_max*p;
            else
                y_min = y_min + y_min*p;
                y_max = y_max + y_min*p;
            end
            ylim(ax,[y_min,y_max]);
            xlim(ax,[0,max(t)]);
            ax.LineWidth = frame_prop;
            title(fig_title,'FontSize',9)
            xlabel('$t$ in $[s]$','FontSize',8)
            ylabel(fig_label,'FontSize',8)
            sp_axes = gca;
            sp_axes.YAxis.Exponent = 0;
            ax.FontSize = 7;
        end
        
        %%
        function fig = plot_staudruck(obj,sol,prob,results_name)
            % Berechnung des Staudrucks
            N = size(sol,1);
            qmax_q_diff = zeros(N,3);
            nof_exc = 0;
            for i = 1:N
                q = (prob.alpha*exp(-prob.beta*sol(i,1))*sol(i,4)^2)/2;
                qmax_q_diff(i,:) = [prob.q_max,q,prob.q_max-q];
                if (q > prob.q_max)
                    nof_exc = nof_exc + 1;
                end
            end
            x = linspace(0,N-1,N);
            % Plot
            fig = figure();
            FigW=16;
            FigH=7;
            plot(x,qmax_q_diff(:,1),'r',x,qmax_q_diff(:,2),'b--');
            xlabel('Diskretisierungspunkte');
            ylabel('$q(v(t),h(t))$');
            set(findall(gcf,'-property','FontSize'),'FontSize',8);
            legend('$q_{\max}$','$q(v(t),h(t))$','Interpreter','latex','Location','southoutside','FontSize',7,'Orientation','horizontal');
%             title(['Staudruck \"Ueberpr\"ufung (',num2str(nof_exc),' \"Uberschreitungen)'],'Interpreter','latex','FontSize',10);
            title(['Staudruck \"Ueberpr\"ufung'],'Interpreter','latex','FontSize',10);
            set(fig,'defaulttextinterpreter','latex',...
                        'PaperUnits','centimeters','PaperSize',[FigW FigH],...
                        'PaperPosition',[0,0,FigW,FigH],'Units','centimeters',...
                        'Position',[0,0,FigW,FigH]);
            set(fig, 'PaperPositionMode', 'auto');
            fig_settings = findobj(fig, 'type', 'axes');
            fig_settings.Units = 'centimeters';
            test=fig_settings.Position;
            fig_settings.Position=[test(1),test(2)-0.35,test(3),test(4)];
            print(fig,'-dpdf','-r600',strcat('./results/',results_name,'_staudruck.pdf'));
        end
    end 
end
