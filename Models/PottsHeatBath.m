classdef PottsHeatBath
    properties
        N
        q
        qs
        beta
        lattice
        J = 1
        E
        fig
        ax
        img
        preCalcExp
    end
    methods
        function obj = PottsHeatBath(N, q, beta)
            obj.N = N;
            obj.q = q;
            obj.beta = beta;
            obj = obj.setup_lattice();
            obj.E = obj.calc_E();
            obj.qs = 1:q;
            obj = obj.setup_figure();
        end
        function obj = set.beta(obj, beta)
            obj.beta = beta;
            obj = obj.update_precalc();
        end
        function obj = update_precalc(obj)
            obj.preCalcExp = exp(obj.beta * obj.J * (0:4));
        end
        function obj = setup_lattice(obj)
            obj.lattice = randi(obj.q, obj.N);
        end
        function E = calc_E(obj)
            b_1 = obj.lattice == circshift(obj.lattice, 1);
            b_2 = obj.lattice == circshift(obj.lattice, 1, 2);
            E = -obj.J * sum(b_1 + b_2, 'all');
        end
        function obj = step(obj)
            y = randi(obj.N);
            x = randi(obj.N);

            y_n1 = mod(y-2, obj.N) + 1;
            x_n1 = mod(x-2, obj.N) + 1;
            y_n2 = mod(y, obj.N) + 1;
            x_n2 = mod(x, obj.N) + 1;

            n_1 = obj.lattice(y_n1, x);
            n_2 = obj.lattice(y, x_n1);
            n_3 = obj.lattice(y_n2, x);
            n_4 = obj.lattice(y, x_n2);

            nbs = [n_1, n_2, n_3, n_4];
            ps = zeros(1, obj.q);
            for s = 1:obj.q
                matchCount = sum(nbs == s);
                ps(s) = obj.preCalcExp(matchCount + 1);
            end
            ps = ps / sum(ps);

            val = rand();
            new_q = find(val <= cumsum(ps));
            obj.lattice(y, x) = new_q(1);
        end
        function obj = setup_figure(obj)
            obj.fig = figure('Visible', 'off');
            obj.ax = axes('Parent', obj.fig);
            obj.img = imagesc(obj.lattice, 'Parent', obj.ax);
            cmap = lines(obj.q);
            colormap(obj.ax, cmap);
            clim(obj.ax, [1 obj.q]);
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