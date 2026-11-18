classdef KawasakiIsing < IsingBase
    properties
        p
    end
    methods
        function obj = KawasakiIsing(n, p, beta)
            obj@IsingBase(n, beta)
            obj.p = p;
            obj = obj.overwrite_lattice();
            obj = obj.calc_E();
            obj = obj.calc_M();
        end
        function obj = overwrite_lattice(obj)
            obj.lattice = ones(obj.n, obj.n);
            n_spins = round(obj.p* obj.n^2);
            locs = randperm(obj.n^2, n_spins);
            obj.lattice(locs) = -1;
        end
        function obj = step(obj)
            x1 = randi(obj.n);
            y1 = randi(obj.n);
            dirs = [1, 0; 0 1; 0 -1; -1 0];
            select_dir = randi(4);
            xdir = dirs(select_dir, 2);
            ydir = dirs(select_dir, 1);
            
            y2 = mod(y1 + ydir - 1, obj.n) + 1;
            x2 = mod(x1 + xdir - 1, obj.n) + 1;

            if obj.lattice(y1, x1) ~= obj.lattice(y2, x2)
                delta_E = 0;
                s1 = obj.lattice(y1, x1);
                s2 = obj.lattice(y2, x2);
                for i = 1:4
                    y2n = mod(y2 + dirs(i, 1) - 1, obj.n) + 1;
                    x2n = mod(x2 + dirs(i, 2) - 1, obj.n) + 1;
                    y1n = mod(y1 + dirs(i, 1) - 1, obj.n) + 1;
                    x1n = mod(x1 + dirs(i, 2) - 1, obj.n) + 1;
                
                    delta_E = delta_E + s2 * obj.lattice(y2n, x2n);
                    delta_E = delta_E + s1 * obj.lattice(y1n, x1n);
                end
                delta_E = delta_E - 2 * s1 * s2;
                delta_E = delta_E * 2 * obj.J;

                if (delta_E <= 0) || (rand() < exp(-obj.beta * delta_E))
                    obj.lattice(y1, x1) = s2;
                    obj.lattice(y2, x2) = s1;
                    obj.E = obj.E + delta_E;
                end
            end
        end
    end
end