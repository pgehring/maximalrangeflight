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
            close("all")
            set(0,'defaulttextinterpreter','latex')
        end
        
        %%
        function fig = plot_fmincon(obj, t, sol, titles, labels, order, frame_prop, line_style)
            obj.nfigures = obj.nfigures +1;
            fig = create_fig(obj, obj.nfigures);
            [r, c] = size(sol);
            for i = [1:c]
                sol_id = order(i);
                obj = create_subplot(obj, t, sol(:,sol_id), [2, 3], i, frame_prop(sol_id), line_style(sol_id));
                title(titles(sol_id))
                xlabel('$t$ in $[s]$')
                ylabel(labels(sol_id))
            end
        end
        
        function fig = create_fig(obj, fig_id)
            fig = figure('WindowState', 'maximized');
            obj.figures(fig_id) = fig;            
        end
        
        function obj = create_subplot(obj, t, y, pos, ax_id, frame_prop, line_style)
            ax = subplot(pos(1), pos(2), ax_id);
            obj.axes(ax_id) = ax;
            if (nargin <= 6)
                plot(t, y, '-k')
            else
                plot(t, y, line_style)
            end
            ax.LineWidth = frame_prop;
        end
    end 
end
