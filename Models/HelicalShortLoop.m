classdef HelicalShortLoop < ShortLoopIceType
    methods
        function obj = HelicalShortLoop(N)
            obj@ShortLoopIceType(N);
        end
        function obj = step(obj)
            visited = false(obj.N^2, 1);
            stack = zeros(obj.N^2, 1);
            ds = zeros(obj.N^2, 1);
            counter = 1;

            % Pick a random lattice site to start the loop
            idx = randi(obj.N^2);
            
            % Add seed lattice site to stack
            stack(counter) = idx;
            
            while ~visited(idx)
                visited(idx) = true;
                [y, x] = ind2sub([obj.N, obj.N], idx);
                % Pick a random outwards direction
                d = obj.lattice(y, x, :);
                d = find(d(:) == 1, 2);
                d = d(randi(2));
                ds(counter) = d;
                
                % Find the neighbour's index
                if d == 1
                    n_idx = mod(idx + obj.N - 1, obj.N^2) + 1;
                elseif d == 2
                    n_idx = mod(idx, obj.N^2) + 1;
                elseif d == 3
                    n_idx = mod(idx - obj.N - 1, obj.N^2) + 1;
                else
                    n_idx = mod(idx - 2, obj.N^2) + 1;
                end
                
                % Add the neighbour to the stack
                counter = counter + 1;
                stack(counter) = n_idx;
                
                % Use the neighbour in the next iteration
                idx = n_idx;
            end
            % Construct loop
            idxs = find(stack == idx);
            start_idx = idxs(1);
            end_idx = idxs(2) - 1;
            
            % Convert to indices for a vectorised update
            latt_idxs = stack(start_idx:end_idx);
            zs = ds(start_idx:end_idx);
            neg_zs = circshift(obj.neg_directions(ds(start_idx:end_idx)), 1)';
            idx_1 = latt_idxs + (zs - 1) * obj.N^2;
            idx_2 = latt_idxs + (neg_zs - 1) * obj.N^2;
            
            % Update edges
            obj.lattice(idx_1) = -obj.lattice(idx_1);
            
            obj.lattice(idx_2) = -obj.lattice(idx_2);
        end
        function result = validate_lattice(obj)
            idxs = 1:obj.N^2;
            cond_1 = obj.lattice(mod(idxs - obj.N - 1, obj.N^2) + 1) == -obj.lattice(idxs + 2*obj.N^2);
            cond_2 = obj.lattice(mod(idxs - 2, obj.N^2) + 1 + obj.N^2) == -obj.lattice(idxs + 3*obj.N^2);
            cond_3 = reshape(sum(obj.lattice, 3) == 0, [1, obj.N^2]);
            result = all(cond_1 & cond_2 & cond_3, 'all');
        end
    end
end