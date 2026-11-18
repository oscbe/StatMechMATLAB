classdef SwendsenWangIsing < IsingBase
    methods
        function obj = SwendsenWangIsing(n, beta)
            obj@IsingBase(n, beta)
        end

        function obj = step(obj)
            p = 1 - exp(-2 * obj.beta * obj.J);
            vertical_bonds = (obj.lattice == circshift(obj.lattice, 1, 1)) .* (rand(obj.n) < p);
            horizontal_bonds = (obj.lattice == circshift(obj.lattice, 1, 2)) .* (rand(obj.n) < p);

            num_clusters = 0;
            checked = zeros(obj.n);
            for x = 1:obj.n
                for y = 1:obj.n
                    if checked(y, x) == 0
                        num_clusters = num_clusters + 1;
                        checked(y, x) = num_clusters;
                        links_to_check = zeros([obj.n * obj.n, 2]);
                        links_to_check(1, 1) = y;
                        links_to_check(1, 2) = x;
                        num_checked = 1;
                        cluster_size = 1;
                        while num_checked <= cluster_size
                            y1 = links_to_check(num_checked, 1);
                            x1 = links_to_check(num_checked, 2);
                            if vertical_bonds(y1, x1) && (checked(mod(y1, obj.n)+1, x1) == 0)
                                checked(mod(y1, obj.n) + 1, x1) = num_clusters;
                                cluster_size = cluster_size + 1;
                                links_to_check(cluster_size, 1) = mod(y1, obj.n) + 1;
                                links_to_check(cluster_size, 2) = x1;
                            end
                            if vertical_bonds(mod(y1-2, obj.n)+1, x1) && (checked(mod(y1-2, obj.n)+1, x1) == 0)
                                checked(mod(y1-2, obj.n)+1, x1) = num_clusters;
                                cluster_size = cluster_size + 1;
                                links_to_check(cluster_size, 1) = mod(y1-2, obj.n) + 1;
                                links_to_check(cluster_size, 2) = x1;
                            end
                            if horizontal_bonds(y1, x1) && (checked(y1, mod(x1, obj.n)+1) == 0)
                                checked(y1, mod(x1, obj.n) + 1) = num_clusters;
                                cluster_size = cluster_size + 1;
                                links_to_check(cluster_size, 1) = y1;
                                links_to_check(cluster_size, 2) = mod(x1, obj.n) + 1;
                            end
                            if horizontal_bonds(y1, mod(x1-2, obj.n)+1) && (checked(y1, mod(x1-2, obj.n)+1) == 0)
                                checked(y1, mod(x1-2, obj.n)+1) = num_clusters;
                                cluster_size = cluster_size + 1;
                                links_to_check(cluster_size, 1) = y1;
                                links_to_check(cluster_size, 2) = mod(x1-2, obj.n) + 1;
                            end
                            num_checked = num_checked + 1;
                        end
                    end
                end
            end
            rand_spins = randi(2, [1, num_clusters])*2-3;
            for i = 1:num_clusters
                obj.lattice(checked == i) = rand_spins(i);
            end
        end
    end
end