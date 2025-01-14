from ea import EA
import sys 

seed = int(sys.argv[1])
periods = [0.0666, 0.05] 

# choose one eval method:

# evaluation_method = "varied"
# evaluation_method = "simult"
evaluation_method = "seq"
# evaluation_method = "tri_obj"

gate_1 = "XOR"
gate_2 = "AND"

source1 = 13
source2 = 41
output = 23

amplitude = 0.005

ea = EA(seed=seed, popsize=100, generations=100, source1=source1, source2=source2, output=output, amplitude=amplitude, periods=periods, evaluation_method=evaluation_method, gate_1=gate_1, gate_2=gate_2)
ea.create_materials()
f = ea.run_generation_one()
ea.hillclimber(f=f)
ea.visualize_best()