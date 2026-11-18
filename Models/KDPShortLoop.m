classdef KDPShortLoop < ShortLoopIceType
    properties
        beta
        E
        eps
    end
    methods
        function obj = KDPShortLoop(N, beta)
            obj@ShortLoopIceType(N)
            obj.beta = beta;
            obj.eps = 1;
            obj.E = obj.calc_E();
        end
        function E = calc_E(obj)
            symmetric = obj.lattice(:, :, 1) == obj.lattice(:, :, 3);
            E = sum(symmetric, 'all');
        end
        function obj = step(obj)
            visited = false(obj.N);
            stack = zeros(obj.N^2, 2);
            ds = zeros(obj.N^2, 1);
            counter = 1;

            % Pick a random lattice site to start the loop
            y = randi(obj.N);
            x = randi(obj.N);
            
            % Add seed lattice site to stack
            stack(counter, 1) = y;
            stack(counter, 2) = x;
            
            while ~visited(y, x)
                visited(y, x) = true;

                % Pick a random outwards direction
                d = obj.lattice(y, x, :);
                d = find(d(:) == 1, 2);
                d = d(randi(2));
                ds(counter) = d;
                
                % Find the neighbour's index
                n_y = obj.wrap(y + obj.directions(d, 1) + 1);
                n_x = obj.wrap(x + obj.directions(d, 2) + 1);
                
                % Add the neighbour to the stack
                counter = counter + 1;
                stack(counter, 1) = n_y;
                stack(counter, 2) = n_x;
                
                % Use the neighbour in the next iteration
                y = n_y;
                x = n_x;
            end
            
            % Construct loop
            idxs = find(all(stack == [y, x], 2));
            start_idx = idxs(1);
            end_idx = idxs(2) - 1;
            
            % Convert to indices for a vectorised update
            ys = stack(start_idx:end_idx, 1);
            xs = stack(start_idx:end_idx, 2);
            zs = ds(start_idx:end_idx);
            neg_zs = circshift(obj.neg_directions(ds(start_idx:end_idx)), 1)';
            idx_1 = sub2ind([obj.N, obj.N, 4], ys, xs, zs);
            idx_2 = sub2ind([obj.N, obj.N, 4], ys, xs, neg_zs);
            
            % Calculate energy change if move is accepted
            energetic = obj.lattice(:, :, 1) == obj.lattice(:, :, 2);
            non_energetic = obj.lattice(:, :, 1) == obj.lattice(:, :, 4);
            v_symmetric = (obj.lattice(:, :, 2) == 1) & (obj.lattice(:, :, 4) == 1);
            h_symmetric = (obj.lattice(:, :, 1) == 1) & (obj.lattice(:, :, 3) == 1);
            same_direction = zs == circshift(zs, -1);
            left_turn = (mod(zs - 2, 4) + 1) == circshift(zs, -1);
            right_turn = (mod(zs, 4) + 1) == circshift(zs, -1);

            latt_idxs = sub2ind([obj.N, obj.N], ys, xs);
            m = sum(energetic(latt_idxs), 'all');
            n = sum(...
                (non_energetic(latt_idxs) & same_direction)...
                | (v_symmetric(latt_idxs) & left_turn)...
                | (h_symmetric(latt_idxs) & right_turn), 'all');

            delta_E = obj.eps * (m - n);
            
            % test_latt = obj.lattice;
            % test_latt(idx_1) = -test_latt(idx_1);
            % test_latt(idx_2) = -test_latt(idx_2);
            % E_1 = sum(obj.lattice(:, :, 1) == obj.lattice(:, :, 2), 'all');
            % E_2 = sum(test_latt(:, :, 1) == test_latt(:, :, 2), 'all');
            % delta_E = -obj.eps*(E_2 - E_1);
            if (delta_E <= 0) || (rand() < exp(-obj.beta * delta_E))
                % Update edges
                obj.lattice(idx_1) = -obj.lattice(idx_1);
                obj.lattice(idx_2) = -obj.lattice(idx_2);
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