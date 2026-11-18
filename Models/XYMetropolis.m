classdef XYMetropolis
    properties
        N
        beta
        lattice
        J = 1
        E
        fig
        ax
        img
    end
    methods
        function obj = XYMetropolis(N, beta)
            obj.N = N;
            obj.beta = beta;
            obj = obj.setup_lattice();
            obj.E = obj.calc_E;
            obj = obj.setup_figure();
        end
        function obj = setup_lattice(obj)
            obj.lattice = rand(obj.N) * 2 * pi;
        end
        function E = calc_E(obj)
            b_1 = cos(obj.lattice - circshift(obj.lattice, 1));
            b_2 = cos(obj.lattice - circshift(obj.lattice, 1, 2));
            E = b_1 + b_2;
        end
        function obj = setup_figure(obj)
            obj.fig = figure('Visible', 'off');
            obj.ax = axes('Parent', obj.fig);
            obj.img = imagesc(obj.lattice, 'Parent', obj.ax);
            clim([0, 2*pi]);
            
            cb = colorbar();
            cbTickVals = [0, pi/2, pi, 3*pi/2, 2*pi];
            cbTickLabels = {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'};
            cb.Ticks = cbTickVals;
            cb.TickLabels = cbTickLabels;
            
            colormap(hsv(256))

            axis(obj.ax, 'equal');
            axis(obj.ax, 'off');
        end

        function obj = plot_lattice(obj)
            set(obj.fig, 'Visible', 'on');
            set(obj.img, 'CData', obj.lattice);
            drawnow;
        end
        function obj = step(obj)
            y = randi(obj.N);
            x = randi(obj.N);

            y_n1 = mod(y - 2, obj.N) + 1;
            y_n2 = mod(y, obj.N) + 1;
            x_n1 = mod(x - 2, obj.N) + 1;
            x_n2 = mod(x, obj.N) + 1;

            ns = [obj.lattice(y, x_n1), obj.lattice(y_n1, x)...
                obj.lattice(y, x_n2), obj.lattice(y_n2, x)];

            new_spin = rand() * 2 * pi;
            delta_E = obj.J * (sum(cos(ns - obj.lattice(y, x)) - cos(ns - new_spin)));

            if (delta_E <= 0) || (rand <= exp(-obj.beta * delta_E))
                obj.lattice(x, y) = new_spin;
            end
        end
    end
end