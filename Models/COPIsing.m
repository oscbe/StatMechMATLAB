classdef COPIsing < IsingBase
    properties
        p
        num_neg
        num_pos
        neg_stack
        pos_stack
    end
    methods
        function obj = COPIsing(n, p, beta)
            obj@IsingBase(n, beta)
            obj.p = p;
            obj = obj.overwrite_lattice();
            obj = obj.calc_E();
            obj = obj.calc_M();
        end
        function obj = overwrite_lattice(obj)
            obj.lattice = ones(obj.n, obj.n);
            obj.num_neg = round(obj.p* obj.n^2);
            obj.num_pos = obj.n^2 - obj.num_neg;
            locs = randperm(obj.n^2, obj.num_neg);
            obj.lattice(locs) = -1;
            [row, col] = find(obj.lattice == -1);
            obj.neg_stack = [row, col];
            [row, col] = find(obj.lattice == 1);
            obj.pos_stack = [row, col];
        end
        function obj = step(obj)
            pos_idx = randi(obj.num_pos);
            neg_idx = randi(obj.num_neg);
            pos_x = obj.pos_stack(pos_idx, 2);
            pos_y = obj.pos_stack(pos_idx, 1);
            neg_x = obj.neg_stack(neg_idx, 2);
            neg_y = obj.neg_stack(neg_idx, 1);
            dirs = [1, 0; 0 1; 0 -1; -1 0];
            
            delta_E = 0;
            for i = 1:4
                ydir = dirs(i, 1);
                xdir = dirs(i, 2);

                py = mod(pos_y + ydir - 1, obj.n) + 1;
                px = mod(pos_x + xdir - 1, obj.n) + 1;
                ny = mod(neg_y + ydir - 1, obj.n) + 1;
                nx = mod(neg_x + xdir - 1, obj.n) + 1;
                
                if ~((py == neg_y) && (px == neg_x))
                    delta_E = delta_E + obj.lattice(py, px);
                end
                if ~((ny == pos_y) && (nx == pos_x))
                    delta_E = delta_E - obj.lattice(ny, nx);
                end
            end
            if (delta_E <= 0) || (rand() < exp(-obj.beta * delta_E))
                obj.lattice(pos_y, pos_x) = -obj.lattice(pos_y, pos_x);
                obj.lattice(neg_y, neg_x) = -obj.lattice(neg_y, neg_x);

                obj.pos_stack(pos_idx, 1) = neg_y;
                obj.pos_stack(pos_idx, 2) = neg_x;
                obj.neg_stack(neg_idx, 1) = pos_y;
                obj.neg_stack(neg_idx, 2) = pos_x;
                
            end
        end
    end
end