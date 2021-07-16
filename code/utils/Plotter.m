classdef Plotter

    properties
        nfigures = 0;
        figures;
        axes;
        n;
        m;
        
    end
    methods
        function obj = Plotter()
            close("all")
            set(0,'defaulttextinterpreter','latex')
        end
        
        function obj = plot_fmincon(obj, t, sol, titles, labels, order)
            obj.nfigures = obj.nfigures +1;
            obj = create_fig(obj, obj.nfigures);
            
            [r, c] = size(sol);
            for i = [1:c]
                sol_id = order(i);
                obj = create_subplot(obj, t, sol(:,sol_id), [2, 3], i);
                title(titles(sol_id))
                ylabel('t in [s]')
                xlabel(labels(sol_id))
            end
        end
        
        function obj = create_fig(obj, fig_id)
            fig = figure('WindowState', 'maximized');
            obj.figures(fig_id) = fig;            

        end
        
        function obj = create_subplot(obj, t, y, pos, ax_id, line_style)
            ax = subplot(pos(1), pos(2), ax_id);
            obj.axes(ax_id) = ax;
            
            if (nargin <= 5)
                plot(t, y, '-k')
            else
                plot(t, y, line_style)
            end
        end
    end 
end
