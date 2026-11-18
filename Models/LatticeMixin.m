classdef LatticeMixin
    properties
        fig
        ax
        img
        lattice
    end
    methods
        function obj = setup_figure(obj)
            obj.fig = figure('Visible', 'off');
            obj.ax = axes('Parent', obj.fig);
            obj.img = imagesc(obj.lattice, 'Parent', obj.ax);
            colormap(obj.ax, [1 1 1; 0 0 0]);
            axis(obj.ax, 'equal');
            axis(obj.ax, 'off');
        end

        function obj = plot_lattice(obj)
            set(obj.fig, 'Visible', 'on');
            set(obj.img, 'CData', obj.lattice);
            drawnow;
        end
    end
end