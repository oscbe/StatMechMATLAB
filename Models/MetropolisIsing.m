classdef MetropolisIsing < IsingBase
    methods
        function obj = MetropolisIsing(n, beta)
            obj@IsingBase(n, beta)
        end

        function obj = step(obj)
            rands = randi(obj.n, 2, 1);
            x = rands(2);
            y = rands(1);
            delta_E = 0;
            directions = [1 0; -1 0; 0 1; 0 -1];
            for i = 1:4
                dx = mod(x + directions(i, 2) - 1, obj.n) + 1;
                dy = mod(y + directions(i, 1) - 1, obj.n) + 1;
                delta_E = delta_E + obj.lattice(dy, dx) * obj.J;
            end
            delta_E = 2 * obj.lattice(y, x) * delta_E;
            if delta_E <= 0
                f = true;
            else
                if rand() < exp(-obj.beta * delta_E)
                    f = true;
                else
                    f = false;
                end
            end
            if f
                obj.lattice(y, x) = -obj.lattice(y, x);
                obj.E = obj.E + delta_E;
                obj.M = obj.M + 2*obj.lattice(y, x);
            end
        end
    end
end