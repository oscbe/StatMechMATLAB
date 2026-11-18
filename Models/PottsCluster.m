classdef PottsCluster < PottsHeatBath
    methods
        function obj = PottsCluster(N, q, beta)
            obj@PottsHeatBath(N, q, beta);
        end
        function obj = step(obj)
            % initialise arrays to track cluster progress
            in_cluster = false(obj.N, obj.N);
            max_stack_size = obj.N * obj.N;
            to_check = zeros(max_stack_size, 2);
            
            % choose a random seed spin
            seed = randi(obj.N, 1, 2);
            
            % define the probability of adding a spin to the cluster
            p = 1 - exp(-obj.beta * obj.J);

            % initialise the stack
            to_check(1, :) = seed;
            in_cluster(seed(1), seed(2)) = true;
            num_checked = 1;
            cluster_size = 1;

            directions = [1 0; -1 0; 0 1; 0 -1];

            while num_checked <= cluster_size
                x = to_check(num_checked, 1);
                y = to_check(num_checked, 2);

                for i = 1:4
                    direction = directions(i, :);
                    adj_x = mod(x + direction(1) - 1, obj.N) + 1;
                    adj_y = mod(y + direction(2) - 1, obj.N) + 1;
                    if ~in_cluster(adj_x, adj_y)
                        if obj.lattice(x, y) == obj.lattice(adj_x, adj_y)
                            if rand(1) < p
                                in_cluster(adj_x, adj_y) = true;
                                cluster_size = cluster_size + 1;
                                to_check(cluster_size, :) = [adj_x, adj_y];
                            end
                        end
                    end
                end
                num_checked = num_checked + 1;
            end
            seed_value = obj.lattice(seed(1), seed(2));
            new_value = mod(seed_value + randi(obj.q - 1) - 1, obj.q) + 1;
            obj.lattice(in_cluster) = new_value;
        end
    end
end