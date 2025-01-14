def get_densities(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    particle_den_map = {}

    x = lines.index("Atoms # sphere\n")
    for line in lines[x + 2:]:
        line_split = line.split()
        particle_den_map[int(line_split[0])] = float(line_split[3])
    return particle_den_map


def read_particle_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    particle_loc_map = {}
    for particle in range(1, 50):
        particle_loc_map[particle] = []

    for line in lines:
        # Skip lines that are not relevant to particle data
        if line.startswith(('ITEM: ATOMS', 'ITEM: TIMESTEP', 'ITEM: NUMBER OF ATOMS', 'ITEM: BOX BOUNDS ')):
            continue

        if len(list(line.split())) < 5:
            continue
        line_list = line.split()
        particle_loc_map[int(line_list[0])].append((float(line_list[1]), float(line_list[2])))

    return particle_loc_map