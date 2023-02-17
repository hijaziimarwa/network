#calculate number of edge in triangle
def launch_multiple_simulations(Nf_list, L_list, nb_list, l0_list,kf_list,ko_list, kc_list,seed_list):
    for Nf in Nf_list:
        for L in L_list:
            for nb in nb_list:
                for l0 in l0_list:
                    for kf in kf_list:
                    	for ko in ko_list:
                    	    for kc in kc_list:
                        	    for seed in seed_list:
                            	    launch_simulation(Nf, L, nb, l0,kf,ko, kc, seed)