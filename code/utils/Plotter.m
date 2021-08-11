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
        end
        
        %%
        function fig = plot_fmincon(obj,t,sol,results_name,titles,labels,order,frame_prop,line_style)
            obj.nfigures = obj.nfigures +1;
            fig = create_fig(obj, obj.nfigures);
            
            fig_title_str = strcat('Versuch: \, ','\verb|',results_name,'|');
            fig_title = sgtitle(fig,fig_title_str,'FontSize',20);
            fig_title.Interpreter = 'latex';
            
            [r,c] = size(sol);
            for i = [1:c]
                sol_id = order(i);
                sp = create_subplot(obj,t,sol(:,sol_id),[2, 3],i,titles(sol_id),labels(sol_id),frame_prop(sol_id),line_style(sol_id));
            end
        end
        
        function fig = plot_state(obj, t, X, titles, labels)
            % function to plot the state vector over time
            
            obj.nfigures = obj.nfigures +1;
            fig = create_fig(obj, obj.nfigures);
            
            
            [r,c] = size(X);
            for i = [1:c]
                sp = create_subplot(obj,t,X(:,i),[c, 1],i,titles(i),labels(i),0.5,"b-");
            end
        end
        
        function fig = create_fig(obj, fig_id)
            fig = figure('WindowState', 'maximized');
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
            ax.LineWidth = frame_prop;
            title(fig_title)
            xlabel('$t$ in $[s]$')
            ylabel(fig_label)
        end
    end 
end
