addpath("/home/oscar/Documents/MATLAB/learning")

a = MetropolisIsing(100, 0.1);
a.plot_lattice();
a.beta = 0.45;
for i = 1:1000
    a = a.evolve(1000);
    a.plot_lattice();
end