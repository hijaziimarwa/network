
from itertools import count
import os
import sys


def create_parameters(Nf, nb, l0, ko, seed):
    suffix = '_Nf%g_nb%d_l0%g_ko%g_seed%d' % (Nf, nb, l0, ko, seed)
    paramfile = 'parameters'+suffix
    dir = 'current'+suffix

    f0 = open('parameters', 'r')
    lines = f0.readlines()
    f0.close()

    fparam = open(paramfile, 'w')

    for line in lines:
        l = line.split()
        if (len(l) > 0):
            current_parameter = l[0]
            if (current_parameter == 'seed'):
                fparam.write(current_parameter+' %d\n' % seed)
            elif (current_parameter == 'num_fibers'):
                fparam.write(current_parameter+' %d\n' % Nf)
            elif (current_parameter == 'num_beads_per_fiber'):
                fparam.write(current_parameter+' %g\n' % nb)
            elif (current_parameter == 'interbead_distance'):
                fparam.write(current_parameter+' %g\n' % l0)
            elif (current_parameter == 'bending_stiffness'):
                fparam.write(current_parameter+' %d\n' % ko)

            else:
                fparam.write(line)

    return paramfile, dir


def launch_simulation(Nf, nb, l0, ko, seed):
    paramfile, dir = create_parameters(Nf, nb, l0, ko, seed)
    # run of a specific node node and exclude node 3 of cluster using sbatch with time 2 days
    os.system('./run_simulation.sh '+paramfile+' '+dir)


def launch_multiple_simulations(Nf_list, nb_list, l0_list, ko_list, seed_list):
    for Nf in Nf_list:
        for nb in nb_list:
            for l0 in l0_list:
                for ko in ko_list:
                    for seed in seed_list:
                        launch_simulation(Nf, nb, l0, ko, seed)


Nf_list = [2450]
nb_list = [19]
l0_list = [5.5]
ko_list = [1]


seed_list = [768]
# print(seed_list)

launch_multiple_simulations(Nf_list, nb_list, l0_list, ko_list, seed_list)

