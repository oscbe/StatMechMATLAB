classdef ShortLoopIceType
    properties
        N
        lattice
        fig
        ax
        img
    end
    properties (Access = protected)
        wrap
        quiver_handles
    end
    properties (Constant, Access = protected)
        directions = [0 1; 1 0; 0 -1; -1 0];
        neg_directions = [3, 4, 1, 2]
    end
    methods
        function obj = ShortLoopIceType(N)
            obj.N = N;
            obj.wrap = [obj.N, 1:obj.N, 1];
            obj = obj.initialise_lattice();
        end
        function obj = initialise_lattice(obj)
            % Represent bonds as pointing outwards (+1) or inwards (-1)
            % for each vertex - note that we are double counting.
            obj.lattice = ones([obj.N, obj.N, 4]);
            obj.lattice(:, :, 3) = -1;
            obj.lattice(:, :, 4) = -1;
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
            
            % Update edges
            obj.lattice(idx_1) = ...
                -obj.lattice(idx_1);
            
            obj.lattice(idx_2) = ...
                -obj.lattice(idx_2);
        end
        
        function result =  validate_lattice(obj)
            cond_1 = obj.lattice(:, :, 2) == -circshift(obj.lattice(:, :, 4), -1, 1);
            cond_2 = obj.lattice(:, :, 1) == -circshift(obj.lattice(:, :, 3), -1, 2);
            cond_3 = sum(obj.lattice, 3) == 0;
            result = all(cond_1 & cond_2 & cond_3, 'all');
        end

        function obj = draw_lattice(obj)
            if isempty(obj.fig) || ~isvalid(obj.fig)
                obj.fig = figure;
                obj.ax = axes('Parent', obj.fig);
                hold(obj.ax, 'on');
                axis(obj.ax, 'equal');
                xlim(obj.ax, [0, obj.N + 1]);
                ylim(obj.ax, [0, obj.N + 1]);
                obj.quiver_handles = cell(1, 4); % Initialize
            end
        
            colours = {'r', 'g', 'b', 'k'};
            dirs = obj.directions;
            scale = 0.3;
            [X, Y] = meshgrid(1:obj.N, 1:obj.N);
        
            for d = 1:4
                L = obj.lattice(:, :, d);
                mask = (L == 1);
                x0 = X(mask);
                y0 = Y(mask);
                dx = dirs(d, 2) * scale * ones(size(x0));
                dy = dirs(d, 1) * scale * ones(size(y0));
        
                if isempty(x0)
                    continue;
                end
        
                % Create or update quiver handle
                if isempty(obj.quiver_handles{d}) || ~isvalid(obj.quiver_handles{d})
                    obj.quiver_handles{d} = quiver(obj.ax, x0, y0, dx, dy, 0, ...
                        'Color', colours{d}, 'MaxHeadSize', 2, 'LineWidth', 1);
                else
                    q = obj.quiver_handles{d};
                    set(q, 'XData', x0, 'YData', y0, 'UData', dx, 'VData', dy);
                end
            end
        
            title(obj.ax, "2D Ice Model Lattice");
        end

        function plot_vertex_types(obj)
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
                1 0 0;    % Type 1 - Red
                0 1 0;    % Type 2 - Green
                0 0 1;    % Type 3 - Blue
                1 1 0;    % Type 4 - Yellow
                0 1 1;    % Type 5 - Cyan
                1 0 1     % Type 6 - Magenta
            ]);
            colorbar(obj.ax, 'Ticks', 1:6, 'TickLabels', {'T1','T2','T3','T4','T5','T6'});
            title(obj.ax, 'Vertex Types');
        end
    end
end