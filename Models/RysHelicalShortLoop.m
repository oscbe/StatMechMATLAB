classdef RysHelicalShortLoop < HelicalShortLoop
    properties
        beta
        eps
    end
    methods
        function obj = RysHelicalShortLoop(N, beta)
            obj@HelicalShortLoop(N)
            obj.beta = beta;
            obj.eps = 1;
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

            symmetric = obj.lattice(:, :, 1) == obj.lattice(:, :, 3);
            m = sum(symmetric(latt_idxs), 'all');
            n = sum((~symmetric(latt_idxs) & (zs ~= circshift(zs, -1))));
            delta_E = obj.eps * (m - n);
            if (delta_E <= 0) || (rand() < exp(-obj.beta * delta_E))
                % Update edges
                obj.lattice(idx_1) = ...
                    -obj.lattice(idx_1);
                
                obj.lattice(idx_2) = ...
                    -obj.lattice(idx_2);
            end
        end
        function plot_symmetric_vertices(obj)
            % Prepare figure
            if isempty(obj.fig) || ~isvalid(obj.fig)
                obj.fig = figure;
                obj.ax = axes;
            else
                cla(obj.ax);
            end
        
            % Get lattice
            L = obj.lattice;
        
            % Flatten direction data
            vertex_types = zeros(obj.N, obj.N);
        
            % Define the six vertex configurations (2-in, 2-out)
            % Each row = [right, up, left, down]
            vertex_lookup = [
                1  1 -1 -1;
                1 -1  1 -1;
                1 -1 -1  1;
               -1  1  1 -1;
               -1  1 -1  1;
               -1 -1  1  1;
            ];
        
            for y = 1:obj.N
                for x = 1:obj.N
                    v = squeeze(L(y, x, :))';
                    for t = 1:6
                        if isequal(v, vertex_lookup(t, :))
                            vertex_types(y, x) = t;
                            break;
                        end
                    end
                end
            end
        
            % Plot
            image(obj.ax, vertex_types);
            axis(obj.ax, 'equal');
            colormap(obj.ax, [
                1 1 1;    % Type 1 - Black
                0 0 0;    % Type 2 - White
                1 1 1;    % Type 3 - White
                1 1 1;    % Type 4 - White
                0 0 0;    % Type 5 - Black
                1 1 1     % Type 6 - White
            ]);
            colorbar(obj.ax, 'Ticks', 1:6, 'TickLabels', {'T1','T2','T3','T4','T5','T6'});
            title(obj.ax, 'Vertex Types');
        end
    end
end