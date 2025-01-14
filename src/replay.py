import pickle 
import random 
import numpy as np 
import os
from ea import EA
from ea import Material 
from ea import Grain

def replay(checkpoint_file, poly=False):

    # Load rng state and afpo state
    with open(checkpoint_file, 'rb') as f:
        ea, rng_state, np_rng_state = pickle.load(f)
    
    os.makedirs(f"../dumps/period{ea.periods[0]}_period{ea.periods[1]}", exist_ok=True)

    for org in ea.population: 
        print(org.fitness1, org.fitness2)
        # if org.fitness2 == 4.6858926596102926: 

            # org.replay_material(poly=poly)

def replay_all_poly(checkpoint_file):

    with open(checkpoint_file, 'rb') as f:
        ea, rng_state, np_rng_state = pickle.load(f)
    os.makedirs(f"../dumps/period{ea.periods[0]}_period{ea.periods[1]}", exist_ok=True)
    ea.visualize_best_poly()

filename = '../checkpoints/simult_eval_run0_period0.0666_period0.05_99gens_AND_AND.p'

replay(filename, poly=True)
# replay_all_poly(filename)