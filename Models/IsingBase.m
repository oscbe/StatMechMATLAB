classdef IsingBase
    properties
        n
        beta
        lattice
        E
        M
        J = 1
        B = 0
        fig
        ax
        img
    end

    methods
        function obj = IsingBase(n, beta)
            obj.n = n;
            obj.beta = beta;
            obj = obj.setup_lattice();
            obj = obj.calc_E();
            obj = obj.calc_M();
            obj = obj.setup_figure();
        end

        function obj = setup_lattice(obj)
            obj.lattice = randi(2, obj.n) * 2 - 3;
        end

        function obj = calc_E(obj)
            val = 0;
            for roll = [-1, 1]
                for dim = [1, 2]
                    val = val + sum(circshift(obj.lattice, roll, dim).*obj.lattice, 'all');
                end
            end
            obj.E = -obj.J*val - obj.B*sum(obj.lattice, 'all');
        end

        function obj = calc_M(obj)
            obj.M = sum(obj.lattice, 'all');
        end

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

        function obj = evolve(obj, steps)
            for i = 1:steps
                obj = obj.step();
            end
        end
    end
    methods (Static)
        function step()
            error("Step not implemented for base class.")
        end
    end
end
