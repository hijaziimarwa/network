from itertools import count
import os
import sys


def create_parameters(Nf, L, nb, l0, kc, seed):

    suffix = '_Nf%d_L%g_nb%d_l0%g_kc%g_seed%d' % (Nf, L, nb, l0, kc, seed)
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
            elif (current_parameter == 'L'):
                fparam.write(current_parameter+' %g\n' % L)
            elif (current_parameter == 'num_beads_per_fiber'):
                fparam.write(current_parameter+' %g\n' % nb)
            elif (current_parameter == 'interbead_distance'):
                fparam.write(current_parameter+' %g\n' % l0)
            elif (current_parameter == 'crosslink_stiffness'):
                fparam.write(current_parameter+' %g\n' % kc)

            else:
                fparam.write(line)

    return paramfile, dir


def launch_simulation(Nf, L, nb, l0, kc, seed):
    paramfile, dir = create_parameters(Nf, L, nb, l0, kc, seed)
    #os.system('sbatch -t "00-05:00:00" -J "%s" -o "logs/%s.log" --export=paramfile=%s,dir=%s ./run_simulation.sh' %(paramfile, paramfile, paramfile, dir))


def launch_multiple_simulations(Nf_list, L_list, nb_list, l0_list, kc_list,seed_list):
    for Nf in Nf_list:
        for L in L_list:
            for nb in nb_list:
                for l0 in l0_list:
                    for kc in kc_list:
                        for seed in seed_list:
                            launch_simulation(Nf, L, nb, l0, kc, seed)


Nf_list = [2000]
L_list = [10]
nb_list = [9]
l0_list = [1]
kc_list = [2]
seed_list = [10]
# print(seed_list)

launch_multiple_simulations(Nf_list, L_list, nb_list, l0_list, kc_list, seed_list)
