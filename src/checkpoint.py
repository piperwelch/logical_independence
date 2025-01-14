import pickle 
import random 
import numpy as np 
import os
from ea import EA
from ea import Material 
from ea import Grain
import sys 

seed = int(sys.argv[1])
gen = int(sys.argv[2])
period1 = 0.0666
period2 = 0.05

def continue_from_checkpoint(checkpoint_file, additional_generations):

    
    # Load rng state and afpo state
    with open(checkpoint_file, 'rb') as f:
        EA, rng_state, np_rng_state = pickle.load(f)

    # Reseed the random number generator
    random.setstate(rng_state)
    np.random.set_state(np_rng_state)
    os.makedirs(f"../dumps/period{period1}_period{period2}", exist_ok=True)
    
    EA.num_generations = EA.generation + additional_generations
    new_fitness_data = np.zeros(shape=(EA.num_generations, EA.popsize, 3))

    print(new_fitness_data.shape)
    new_fitness_data[:EA.generation, :, :] = EA.fitness_data[:EA.generation,:,:]

    EA.fitness_data = new_fitness_data
    EA.hillclimber(continue_from_checkpoint=True)

eval_method = 'tri_obj'

checkpoint_filename = f'../checkpoints/{eval_method}_eval_run{seed}_period{period1}_period{period2}_{gen}gens_XOR_AND.p'
# CGMMs/logical_independence/high_fidelity_physics/multiscale_varied/checkpoints/seq_eval_run0_period0.0666_period0.05_52gens_AND_AND.p
continue_from_checkpoint(checkpoint_filename, additional_generations=23)