import random
import numpy as np 
import constants as c 
import os 
import re 
from datetime import datetime
import copy 
import time
import pickle 
import random 
import itertools
from glob import glob 
import subprocess 
from itertools import product
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt


class Grain:
    def __init__(self, x, y, z, density, grain_id, is_source=False, is_end=False):
        self.x = x
        self.y = y
        self.z = z
        self.density = density
        self.is_source = is_source
        self.is_end = is_end
        self.grain_id = grain_id


    def mutate(self, mutation_rate=0.1, max_mutation_amount=100):
        # Mutate the density of the grain with a certain probability
        if random.random() < mutation_rate:
            mutation_amount = random.uniform(-max_mutation_amount, max_mutation_amount)
            self.density += mutation_amount
        if self.density < 500: #kg/m3
            self.density = 500


class Material:
    def __init__(self, id, periods, source1, source2, amplitude, output, eval_method, gate_1, gate_2, seed = 0):
        self.id = id
        self.seed = seed
        self.evaluation_method = eval_method
        if eval_method == "tri_obj": self.conditions = ['0101', '1010', '1111']
        if eval_method == "simult": self.conditions = ['0101', '1010', '1111']
        if eval_method == "seq": self.conditions = ['0001', '0010', '1000', '0100', '1100', '0011']
        self.gate_1 = gate_1 
        self.gate_2 = gate_2 

        self.amplitude = amplitude
        self.periods = periods
        self.source1 = source1
        self.source2 = source2
        self.output = output
        self.grain_diameter = 0.1
        self.make_grains()


    def gen_position_fit_dict(self):
        '''generates input case fitness dictionary'''
        self.period_condition_fit_dict = {
            condition: {period: -1 for period in self.periods}
            for condition in self.conditions
        }
        

    def read_sim(self, condition):
        '''parses LAMMPS dump file'''
        grain_ids = {self.source1, self.source2, self.output}
        positions = {grain_id: [] for grain_id in grain_ids}

        with open(f"../dumps/period{self.periods[0]}_period{self.periods[1]}/dump.{self.evaluation_method}_eval_seed{self.seed}_id{self.id}_{condition}", 'r') as file:
            data=file.read()
            # Split data string into sections based on "ITEM: ATOMS" lines
            sections = data.strip().split("ITEM: ATOMS id x y\n")[1:]

            for section in sections:
                lines = section.strip().split('\n')
                for line in lines:
                    parts = line.strip().split()
                    if len(parts) == 3:
                        try:
                            grain_id = int(parts[0])
                            if grain_id in grain_ids:
                                x, y = map(float, parts[1:3])

                                positions[grain_id].append(x)
                        except:
                            continue
        return positions
    

    def analyze_sim(self, condition, xs, driving_frequency1=None, driving_frequency2=None):
        '''Parse time domain data into frequency domain data'''
        xs = np.array(xs) - np.mean(xs)
        time_data = np.linspace(0, 1, len(xs))
        freq_signal = np.abs(np.fft.rfft(xs))
        frequencies = np.fft.rfftfreq(len(time_data), 0.005)
        
        if driving_frequency1 != None:
            peak_at_f1 = freq_signal[np.abs(frequencies - driving_frequency1).argmin()]
            self.period_condition_fit_dict[condition][self.periods[0]] = peak_at_f1

        if driving_frequency2 != None:
            peak_at_f2 = freq_signal[np.abs(frequencies - driving_frequency2).argmin()]
            self.period_condition_fit_dict[condition][self.periods[1]] = peak_at_f2
    


    def make_grains(self): 
        '''Instantiate material's grains'''
        r = self.grain_diameter/2

        grain_id = 1
        self.grains = []
        for i in range(7):
            for j in range(7):
                for k in range(1):
                    x = (r * (2*i + ((j + k) % 2))) * 0.9
                    y = (r * (np.sqrt(3) * (j + (1/3) * (k % 2)))) *0.9
                    z = r * ((2 * np.sqrt(6) / 3) * k) 

                    grain_density = random.randint(1000, 2500)
                    is_input, is_output  = False, False
                    if grain_id in [self.source1, self.source2]:
                        is_input = True
                    if grain_id == self.output: 
                        is_output = True
     
                    grain = Grain(grain_id=grain_id, density=grain_density, x=x, y=y, z=z, is_source=is_input, is_end = is_output)
                    self.grains.append(grain)
                    grain_id+=1


    def write_data(self):
        '''Write data file for LAMMPS to parse'''

        # f = open(f"../results/{self.evaluation_method}_eval_data_id{self.id}_seed{self.seed}_period{self.periods[0]}_period{self.periods[1]}", "w")
        f = open(f"../data_ins/{self.evaluation_method}_eval_data_id{self.id}_seed{self.seed}_period{self.periods[0]}_period{self.periods[1]}", "w")

        f.write(f"\n{49} atoms\n")
        f.write(f"{2} atom types\n\n")
        f.write(f"{-self.grain_diameter} 0.8 xlo xhi\n")
        f.write(f"{-self.grain_diameter} 0.8 ylo yhi\n")
        f.write(f"{-self.grain_diameter/2} {3*(self.grain_diameter/2)} zlo zhi\n")
        f.write(f"0 0 0 xy xz yz\n\n")
        f.write(f"Atoms # sphere\n\n")

        grain_count = 1
        for grain in self.grains:
            #                id    type radius density           x        y      z
            f.write(f"{grain_count} 1     0.1 {grain.density} {grain.x} {grain.y} {grain.z} \n")
            grain_count+=1 

        f.close()


    def get_fitness(self, evaluation_method):
        '''Get fitness using on appropriate method'''
        if evaluation_method == "varied":
            self.get_fitness_varied()
        if evaluation_method == "tri_obj":
            self.get_fitness_tri_object()
        if evaluation_method == "simult":
            self.get_fitness_simult()
        if evaluation_method == "seq":
            self.get_fitness_seq()


    def get_fitness_tri_object(self):
        fitness1_11 = 0
        fitness1_01 = 0
        fitness1_10 = 0
        fitness2_11 = 0
        fitness2_01 = 0
        fitness2_10 = 0

        for case,behavior in self.period_condition_fit_dict.items():
            if case == '1010':
                fitness1_10+=behavior[self.periods[0]]
                fitness2_10+=behavior[self.periods[1]]
            if case == '0101':
                fitness1_01+=behavior[self.periods[0]]
                fitness2_01+=behavior[self.periods[1]]
            if case == '1111':
                fitness1_11+=behavior[self.periods[0]]
                fitness2_11+=behavior[self.periods[1]]
            
        if self.gate_1 == "AND":
            self.fitness1 = fitness1_11/((fitness1_01 + fitness1_10)/2)
        if self.gate_1 == "XOR":   
            self.fitness1 = ((fitness1_01 + fitness1_10)/2)/fitness1_11

        if self.gate_2 == "AND":
            self.fitness2 = fitness2_11/((fitness2_01 + fitness2_10)/2)
        if self.gate_2 == "XOR":
            self.fitness2 = ((fitness2_01 + fitness2_10)/2)/fitness2_11
        
        if self.gate_1 == "AND" and self.gate_2 == "AND":
            signal_difference = np.abs(fitness2_11 - fitness1_11)
        elif self.gate_1 == "XOR":
            signal_difference = np.abs(fitness2_11 - (fitness1_01+fitness1_10)/2)
        else:
            signal_difference = np.abs(fitness1_11 - (fitness2_01+fitness2_10)/2)

        weight = 3 
        self.fitness1*=1/(1+(signal_difference*weight))
        self.fitness2*=1/(1+(signal_difference*weight))


    def get_fitness_simult(self):
        fitness1_11 = 0
        fitness1_01 = 0
        fitness1_10 = 0
        fitness2_11 = 0
        fitness2_01 = 0
        fitness2_10 = 0

        for case,behavior in self.period_condition_fit_dict.items():
            if case == '1010':
                fitness1_10+=behavior[self.periods[0]]
                fitness2_10+=behavior[self.periods[1]]
            if case == '0101':
                fitness1_01+=behavior[self.periods[0]]
                fitness2_01+=behavior[self.periods[1]]
            if case == '1111':
                fitness1_11+=behavior[self.periods[0]]
                fitness2_11+=behavior[self.periods[1]]
            
        if self.gate_1 == "AND":
            self.fitness1 = fitness1_11/((fitness1_01 + fitness1_10)/2)
        if self.gate_1 == "XOR":   
            self.fitness1 = ((fitness1_01 + fitness1_10)/2)/fitness1_11

        if self.gate_2 == "AND":
            self.fitness2 = fitness2_11/((fitness2_01 + fitness2_10)/2)
        if self.gate_2 == "XOR":
            self.fitness2 = ((fitness2_01 + fitness2_10)/2)/fitness2_11
        
        if self.gate_1 == "AND" and self.gate_2 == "AND":
            signal_difference = np.abs(fitness2_11 - fitness1_11)
        elif self.gate_1 == "XOR":
            signal_difference = np.abs(fitness2_11 - (fitness1_01+fitness1_10)/2)
        else:
            signal_difference = np.abs(fitness1_11 - (fitness2_01+fitness2_10)/2)

        weight = 3 
        self.fitness1*=1/(1+(signal_difference*weight))
        self.fitness2*=1/(1+(signal_difference*weight))


    def get_threshold(self,independence_results): 
        highest_zero = -1
        lowest_one = 1000
        for input_case_pair, behavior in independence_results.items(): 
            input_f1 = input_case_pair[:2]
            input_f2 = input_case_pair[2:]
            behavior_f1 = behavior[0]
            behavior_f2 = behavior[1]

            if input_f1 == '11' and behavior_f1 < lowest_one: 
                lowest_one = behavior_f1
            if input_f2 == '11' and behavior_f2 < lowest_one: 
                lowest_one = behavior_f2
            if input_f1 != '11' and behavior_f1 > highest_zero: 
                highest_zero = behavior_f1
            if input_f2 != '11' and behavior_f2 > highest_zero: 
                highest_zero = behavior_f2

        return (highest_zero+lowest_one)/2


    def analyze_sim_independence(self, condition, xs, driving_frequency1, driving_frequency2):
        #change to displacement 
        xs = np.array(xs) - np.mean(xs)
        time_data = np.linspace(0, 1, len(xs))
        freq_signal = np.abs(np.fft.rfft(xs))
        frequencies = np.fft.rfftfreq(len(time_data), 1)
        
        peak_at_f1 = freq_signal[np.abs(frequencies - driving_frequency1).argmin()]
        peak_at_f2 = freq_signal[np.abs(frequencies - driving_frequency2).argmin()]

        return peak_at_f1, peak_at_f2

    def get_independence(self, independence_cases, driving_frequency1=15, driving_frequency2=20):
        independence_results = {}

        #get peaks from the sims 
        for input_case_pair in independence_cases: 
            positions = self.read_sim(input_case_pair)[self.output]
            independence_results[input_case_pair] = self.analyze_sim_independence(input_case_pair, positions, driving_frequency1, driving_frequency2)
        #get threshold 
        threshold = self.get_threshold(independence_results)

        #gate_1 is XOR 
        #gate_2 is AND
        independence = 0 
        for input_case_pair, behavior in independence_results.items():
            input_f1 = input_case_pair[:2]
            input_f2 = input_case_pair[2:]
            behavior_f1 = behavior[0]
            behavior_f2 = behavior[1]
            correct_f1, correct_f2 = False, False
            if input_f1 == '01' and behavior_f1 >= threshold: 
                correct_f1 = True
            elif input_f1 == '10' and behavior_f1 >= threshold: 
                correct_f1 = True

            elif input_f1 != '01' and behavior_f1 <= threshold: 
                correct_f1 = True
            elif input_f1 != '10'  and behavior_f1 <= threshold: 
                correct_f1 = True
            else: 
                pass
            
            if input_f2 == '11' and behavior_f2 >= threshold: 
                correct_f2 = True
            if input_f2 != '11' and behavior_f2 <= threshold: 
                correct_f2 = True

            independence += (correct_f1 and correct_f2)
        self.independence = independence


    def get_fitness_seq(self):
        fitness1_11 = 0
        fitness1_01 = 0
        fitness1_10 = 0
        fitness2_11 = 0
        fitness2_01 = 0
        fitness2_10 = 0

        for case,behavior in self.period_condition_fit_dict.items():
            if case == '0001':
                fitness2_01+=behavior[self.periods[1]]
            if case == '0010':
                fitness2_10+=behavior[self.periods[1]]
            if case == '0011':
                fitness2_11+=behavior[self.periods[1]]
            if case == '0100':
                fitness1_01+=behavior[self.periods[0]]
            if case == '1000':
                fitness1_10+=behavior[self.periods[0]]
            if case == '1100':
                fitness1_11+=behavior[self.periods[0]] 


        if self.gate_1 == "AND":
            self.fitness1 = fitness1_11/((fitness1_01 + fitness1_10)/2)
        if self.gate_1 == "XOR":   
            self.fitness1 = ((fitness1_01 + fitness1_10)/2)/fitness1_11

        if self.gate_2 == "AND":
            self.fitness2 = fitness2_11/((fitness2_01 + fitness2_10)/2)
        if self.gate_2 == "XOR":
            self.fitness2 = ((fitness2_01 + fitness2_10)/2)/fitness2_11

        if self.gate_1 == "AND" and self.gate_2 == "AND":
            signal_difference = np.abs(fitness2_11 - fitness1_11)
        elif self.gate_1 == "XOR":
            signal_difference = np.abs(fitness2_11 - (fitness1_01+fitness1_10)/2)
        else:
            signal_difference = np.abs(fitness1_11 - (fitness2_01+fitness2_10)/2)

        self.fitness1*=1/(1+(signal_difference))
        self.fitness2*=1/(1+(signal_difference))


    def get_fitness_varied(self):
    
        fitness1_11 = 0
        fitness1_01 = 0
        fitness1_10 = 0
        fitness1_count_11 = 0 
        fitness1_count_01 = 0 
        fitness1_count_10 = 0 

        fitness2_11 = 0
        fitness2_01 = 0
        fitness2_10 = 0
        fitness2_count_11 = 0 
        fitness2_count_01 = 0 
        fitness2_count_10 = 0 


        for case,behavior in self.period_condition_fit_dict.items():
            if case[:2] == "11":
                fitness1_11+=behavior[self.periods[0]]
                fitness1_count_11+=1
            if case[:2] == "01":
                fitness1_01+=behavior[self.periods[0]]
                fitness1_count_01+=1 
            if case[:2] == "10":
                fitness1_10+=behavior[self.periods[0]]
                fitness1_count_10+=1 

            if case[2:] == "11":
                fitness2_11+=behavior[self.periods[1]]
                fitness2_count_11+=1
            if case[2:] == "01":
                fitness2_01+=behavior[self.periods[1]]
                fitness2_count_01+=1 
            if case[2:] == "10":
                fitness2_10+=behavior[self.periods[1]]
                fitness2_count_10+=1 
        
        #make for safe division
        peak_term = False
        needs2_AND_div,needs1_AND_div,needs2_XOR_div,needs1_XOR_div = False,False,False,False
        if fitness2_11 != 0 and fitness1_11 != 0:
            peak_term = True; signal_difference = np.abs(fitness2_11 - fitness1_11)


        if self.gate_1 == "AND" and self.gate_2 == "AND" and fitness2_11 != 0 and fitness1_11 != 0:
            signal_difference = np.abs(fitness2_11 - fitness1_11)
            peak_term = True
        if self.gate_1 == "XOR" and (fitness1_01 != 0 and fitness1_10 != 0) and fitness2_11 !=0:
            signal_difference = np.abs(fitness2_11 - (fitness1_01+fitness1_10)/2)
            peak_term = True
        if self.gate_2 == 'XOR' and (fitness2_01 != 0 and fitness2_10 != 0) and fitness1_11 !=0:
            signal_difference = np.abs(fitness1_11 - (fitness2_01+fitness2_10)/2)
            peak_term = True

        if fitness1_count_11 == 0: 
            fitness1_count_11 = 1
            fitness1_11 = 1

        if fitness2_count_11 == 0: 
            fitness2_count_11 = 1
            fitness2_11 = 1

        if fitness1_count_01 == 0 and fitness1_count_10 == 0:
            needs1_AND_div = True 
            needs1_XOR_div = True 

        if fitness2_count_01 == 0 and fitness2_count_10 == 0:
            needs2_AND_div = True 
            needs2_XOR_div = True 

   
        if self.gate_1 == "AND": 
            if needs1_AND_div: 
                self.fitness1 = fitness1_11/fitness1_count_11
            else:
                self.fitness1 = (fitness1_11/fitness1_count_11)/((fitness1_01 + fitness1_10)/(fitness1_count_01+fitness1_count_10))
        elif self.gate_1 == "XOR":
            if needs1_XOR_div:
                self.fitness1 = 1/(fitness1_11/fitness1_count_11)
            else: 
                self.fitness1 = ((fitness1_01 + fitness1_10)/(fitness1_count_01+fitness1_count_10))/(fitness1_11/fitness1_count_11)

        if self.gate_2 == "AND":
            if needs2_AND_div: 
                self.fitness2 = fitness2_11/fitness2_count_11
            else:
                self.fitness2 = (fitness2_11/fitness2_count_11)/((fitness2_01 + fitness2_10)/(fitness2_count_01+fitness2_count_10))
        elif self.gate_2 == "XOR": 
            if needs2_XOR_div: 
                self.fitness2 = 1/(fitness2_11/fitness2_count_11)
            else:
                self.fitness2 = ((fitness2_01 + fitness2_10)/(fitness2_count_01+fitness2_count_10))/(fitness2_11/fitness2_count_11)
       
        if peak_term: 
            self.fitness1*= 1/(1+signal_difference)
            self.fitness2*= 1/(1+signal_difference)


    def mutate(self):
        ''' mutates all grains in a material 
        '''
        for grain in self.grains:
            grain.mutate()


    def replay_material(self, replay=True, poly=True):
        '''Replays materials for all 15 input cases'''

        os.makedirs(f"../dumps/period{self.periods[0]}_period{self.periods[1]}", exist_ok=True)

        self.write_data()
        current_datetime = datetime.now()
        if replay:
            data_file = f"../results/{self.evaluation_method}_eval_data_id{self.id}_seed{self.seed}_period{self.periods[0]}_period{self.periods[1]}"

        else:
            data_file = f"../data_ins/{self.evaluation_method}_eval_data_id{self.id}_seed{self.seed}_period{self.periods[0]}_period{self.periods[1]}"
        for condition in ["0001", "0010", "0011", "0100", "1000", "1100", "0110", "1001", "1111", "0111", "1011", "1101", "1110", "0101", "1010"]: 
           os.system(f"./lmp_serial -in in.mat_simulator_varied -var data_file {data_file} \
            -var source_node1 {self.source1} -var source_node2 {self.source2} -var output_node {self.output} -var id {self.id} \
            -var visualize True -var seed {self.seed} -var amplitude {self.amplitude} -var period1 {self.periods[0]} -var period2 {self.periods[1]} \
            -var method {self.evaluation_method} -var case {condition} -var dump_folder ../dumps/period{0.0666}_period{0.05}/ >../outputs/{11_14_2024}")
           


class EA:
    def __init__(self, seed, popsize, generations, periods, source1, source2, output, amplitude, evaluation_method, gate_1, gate_2, continue_from_checkpoint=False):
        self.seed = seed 
        self.popsize = popsize 
        self.num_generations = generations 
        self.population = []
        self.next_avaliable_id = 0 
        self.periods = periods
        self.source1 = source1
        self.source2 = source2
        self.output = output 
        self.amplitude = amplitude
        self.generation = 0 
        self.evaluation_method = evaluation_method
        self.gate_1 = gate_1
        self.gate_2 = gate_2 
        os.makedirs(f"../dumps/period{self.periods[0]}_period{self.periods[1]}/", exist_ok=True)
        np.random.seed(seed)
        random.seed(seed)
        if self.evaluation_method == "tri_obj":
            self.fitness_data = np.zeros(shape=(self.num_generations+1, self.popsize, 3))
        else:
            self.fitness_data = np.zeros(shape=(self.num_generations+1, self.popsize, 2))
    

    def generate_independence_metrics(self): 

        bit_strings = [''.join(bits) for bits in itertools.product('01', repeat=4)]

        # Define a subset of bit strings to remove
        subset_to_remove = {'0000', '1111', '0101', '1010'}

        # Remove the specified subset from the list of bit strings
        filtered_bit_strings = [bit_string for bit_string in bit_strings if bit_string not in subset_to_remove]
        selected_bit_strings = []
        while True:
            selected_bit_strings = random.sample(filtered_bit_strings, 3)
            if any(bit_string[:2] == '11' or bit_string[-2:] == '11'for bit_string in selected_bit_strings):
                break
        self.indepence_cases = selected_bit_strings


    def create_materials(self):
        for _ in range(self.popsize*2):
            self.population.append(Material(id=self.next_avaliable_id, periods=self.periods, \
            source1=self.source1, source2=self.source2, amplitude=self.amplitude, \
            output=self.output, seed=self.seed, eval_method=self.evaluation_method, gate_1 = self.gate_1, gate_2 = self.gate_2))
            self.next_avaliable_id+=1 


    def run_batch(self, organisms): 
        data_files = ""
        all_cases = ""
        for organism in organisms: 
            data_file = f"../data_ins/{self.evaluation_method}_eval_data_id{organism.id}_seed{self.seed}_period{organism.periods[0]}_period{organism.periods[1]} "
            cases = ",".join(organism.conditions)
            data_files += data_file
            all_cases += cases + " "
        if self.evaluation_method == "tri_obj":
            cases +=  f",{self.indepence_cases[0]}" + f",{self.indepence_cases[1]}" + f",{self.indepence_cases[2]}" + " "
            all_cases = cases*10

        current_datetime = datetime.now()
        os.system(f"sbatch batch10_run_condition.sh {data_files} {self.source1} {self.source2} {self.output} False {self.amplitude} {self.evaluation_method} {current_datetime} {all_cases}")


    def run_batch_direct(self, organisms):
        current_datetime = datetime.now()

        for organism in organisms:
            data_file = f"../data_ins/{self.evaluation_method}_eval_data_id{organism.id}_seed{self.seed}_period{organism.periods[0]}_period{organism.periods[1]} "
            
            if self.evaluation_method == "tri_obj":
                curr_conditions = organism.conditions + self.indepence_cases
            else:
                curr_conditions = organism.conditions 
            for case in curr_conditions:
                os.system(f"./lmp_serial -in in.mat_simulator_varied -var data_file {data_file} \
                -var source_node1 {self.source1} -var source_node2 {self.source2} \
                -var output_node {self.output} -var id {organism.id} -var visualize {False} \
                -var seed {self.seed} -var amplitude {self.amplitude} -var period1 {organism.periods[0]} \
                -var period2 {organism.periods[1]}  -var method {self.evaluation_method} -var case {case} \
                -var dump_folder ../dumps/period{organism.periods[0]}_period{organism.periods[1]}/ \
                >../outputs/output &")


    def generate_random_case(self):
        cases = ["00", "10", "01", "11"]
        chosen_cases = []
        for _ in range(3):
            potential_case = ''.join(random.choices(cases, k=2))
            while potential_case == '0000' or potential_case in chosen_cases:
                potential_case = ''.join(random.choices(cases, k=2))
            chosen_cases.append(potential_case)
        return chosen_cases


    def run_generation_one(self):
        if self.evaluation_method == "tri_obj": 
            self.generate_independence_metrics()

        for material in self.population:
            material.write_data()
        if self.evaluation_method == "varied":
            curr_gen_cases = self.generate_random_case()
            for material in self.population: 
                material.conditions = curr_gen_cases
        for material in self.population: material.gen_position_fit_dict()

        for i in range(0, len(self.population), 10):
            run_batch = self.population[i:i+10]  # Get the segment of length 10
            self.run_batch_direct(run_batch)  # Pass the segment to another functi
  
        self.wait_for_sims_to_finish()

        for material in self.population:
            for condition in material.conditions:
                positions = material.read_sim(condition)
                driving_freq1=15
                driving_freq2=20
                material.analyze_sim(condition, positions[self.output], driving_freq1, driving_freq2)
            
            material.get_fitness(self.evaluation_method)
            if self.evaluation_method == "tri_obj":
                material.get_independence(self.indepence_cases, driving_freq1, driving_freq2)

        self.generation+=1

        os.system(f"rm -rf ../data_ins/*{self.evaluation_method}_eval*seed{self.seed}_period{self.periods[0]}_period{self.periods[1]}")
        os.system(f"rm -rf ../dumps/period{self.periods[0]}_period{self.periods[1]}/*{self.evaluation_method}_eval_seed{self.seed}*")
        os.system(f"rm outputs/*")


    def wait_for_sims_to_finish(self):

        all_sims_started = False
        if self.evaluation_method == "seq" or self.evaluation_method == "tri_obj":
            n_sims = 6*len(glob(f"../data_ins/*seed{self.seed}*period{self.periods[0]}_period{self.periods[1]}*"))
        else:
            n_sims = 3*len(glob(f"../data_ins/*seed{self.seed}*period{self.periods[0]}_period{self.periods[1]}*"))

        while not all_sims_started:
            total_finished_sims = len(glob(f'../dumps/period{self.periods[0]}_period{self.periods[1]}/*seed{self.seed}*')) 
            if total_finished_sims == n_sims:
                all_sims_started = True
            else: #must update to check appropriate period 
                time.sleep(10) # check in increments of 1 seconds
                print("Seed:", self.seed, "Generation:", self.generation, "Generation Progress:",len(glob(f'../dumps/period{self.periods[0]}_period{self.periods[1]}/*seed{self.seed}*')), "complete of", n_sims, flush=True)
        finished_sims = []
        while len(finished_sims) != n_sims:
            for out_file in glob(f'../dumps/period{self.periods[0]}_period{self.periods[1]}/*seed{self.seed}*'):
                
                f = open(out_file, 'r')
                lines = f.readlines()

                sim_name = out_file.split(f'/dumps/period{self.periods[0]}_period{self.periods[1]}/')[1]
                if sim_name not in finished_sims:

                    if len(lines) > 4000:
                        finished_sims.append(sim_name)

            time.sleep(1) # check in increments of 1 seconds   
            # print("lines 253", n_sims, len(finished_sims))      
    

    def save_checkpoint(self, j):

        filename = '../checkpoints/{}_eval_run{}_period{}_period{}_{}gens_{}_{}.p'.format(self.evaluation_method, self.seed, self.periods[0], self.periods[1], j, self.gate_1, self.gate_2)

        rng_state = random.getstate()
        np_rng_state = np.random.get_state()

        with open(filename, 'wb') as f:
            pickle.dump([self, rng_state, np_rng_state], f)


    def survivor_selection(self):
        '''Tournament selection'''
        tri_obj = True if self.evaluation_method == "tri_obj" else False
        while len(self.population) > self.popsize:

            # Choose two different individuals from the population
            ind1 = np.random.randint(len(self.population))
            ind2 = np.random.randint(len(self.population))
            while ind1 == ind2:
                ind2 = np.random.randint(len(self.population))

            if self.dominates(ind1, ind2, tri_obj=tri_obj):  # ind1 dominates
                # remove ind2 from population and shift following individuals up in list
                for i in range(ind2, len(self.population)-1):
                    self.population[i] = self.population[i+1]
                self.population.pop() # remove last element from list (because it was shifted up)

            elif self.dominates(ind2, ind1, tri_obj=tri_obj):  # ind2 dominates

                # remove ind1 from population and shift following individuals up in list
                for i in range(ind1, len(self.population)-1):
                    self.population[i] = self.population[i+1]
                self.population.pop() # remove last element from list (because it was shifted up)

        assert len(self.population) == self.popsize


    def dominates(self, ind1, ind2, tri_obj=False):
        '''Returns true is individual 1 dominates individual 2'''
        if tri_obj == False: 
            if self.population[ind1].fitness1 >= self.population[ind2].fitness1 and self.population[ind1].fitness2 >= self.population[ind2].fitness2:
                return True
            else:
                return False
        else: 
            if self.population[ind1].fitness1 >= self.population[ind2].fitness1 and self.population[ind1].fitness2 >= self.population[ind2].fitness2 and self.population[ind1].independence >= self.population[ind2].independence:
                return True
            else:
                return False


    def hillclimber(self, continue_from_checkpoint=False, f = None):
        '''Hillclimber'''

        for j in range(self.generation, self.num_generations):
            self.survivor_selection()

            for i, mat in enumerate(self.population):   
                self.fitness_data[self.generation,i,0] = mat.fitness1
                self.fitness_data[self.generation,i,1] = mat.fitness2

                if self.evaluation_method == "tri_obj":
                    self.fitness_data[self.generation,i,1] = mat.independence


            if j % 1 == 0: 
                self.save_checkpoint(j)
                print(j)

            new_orgs = []
            num_children_per_parent = 2 
            #tournament winners breed 
            for parent in self.population:
                for m in range(num_children_per_parent): #for the num children per parent
                    child = copy.deepcopy(parent)
                    child.id = self.next_avaliable_id
                    self.next_avaliable_id+=1
                    child.mutate()
                    child.write_data()
                    new_orgs.append(child)
                if self.evaluation_method == "varied":
                    parent.write_data()
            self.population.extend(new_orgs)

            if self.evaluation_method == "varied":
                curr_gen_cases = self.generate_random_case()
                for material in self.population: 
                    material.conditions = curr_gen_cases
                    material.gen_position_fit_dict()
                    eval_pop = self.population 
            else: 
                eval_pop = new_orgs
            for i in range(0, len(eval_pop), 10):
                run_batch = eval_pop[i:i+10]  # Get the segment of length 10
                self.run_batch_direct(run_batch)  # Pass the segment to another functi

            self.wait_for_sims_to_finish()

            for material in eval_pop:
                for condition in material.conditions:
                    positions = material.read_sim(condition)
                    driving_freq1=15
                    driving_freq2=20
                    material.analyze_sim(condition, positions[self.output], driving_freq1, driving_freq2)
                material.get_fitness(self.evaluation_method)
                if self.evaluation_method == "tri_obj":
                    material.get_independence(self.indepence_cases, driving_freq1, driving_freq2)

            self.generation = j

            os.system(f"rm -rf ../data_ins/*{self.evaluation_method}_eval*_seed{self.seed}_period{self.periods[0]}_period{self.periods[1]}")
            os.system(f"rm -rf ../dumps/period{self.periods[0]}_period{self.periods[1]}/*{self.evaluation_method}_eval_seed{self.seed}*")
            os.system(f"rm outputs/*")


    def visualize_best(self, poly=False):
        top = max(self.population, key=lambda x: x.fitness)
        top.replay_material(poly=poly)


    def visualize_best_poly(self):
        top = max(self.population, key=lambda x: x.fitness)
        data_file = f"../data_ins/data_id{top.id}_seed{self.seed}_period{top.periods[0]}_period{top.periods[1]}"
        os.system(f"sbatch check_indp.sh {data_file} {top.source1} {top.source2} {top.output} {top.amplitude}")