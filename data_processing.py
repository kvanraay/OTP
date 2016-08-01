import world
import csv
from cell import Bacteria, Cell, BACTERIA_TYPES
from collections import Counter

def output_counts_to_csv(file_name, counters):
    """
    writes data from simulation to csv file format
    """
    with open(file_name, 'w') as csvfile:
        columns = ['Update', 'Strain', 'Count']

        writer = csv.DictWriter(csvfile, fieldnames=columns)

        writer.writeheader()

        for update, counter in enumerate(counters):
            for bac_type in sorted(list(BACTERIA_TYPES)):
                if bac_type not in counter:
                    counter[bac_type] = 0

                writer.writerow({"Update": update, "Strain":bac_type, "Count":counter[bac_type]})


        #writer.writerow({'first_name': 'Baked', 'last_name': 'Beans'})
        #writer.writerow({'first_name': 'Lovely', 'last_name': 'Spam'})
        #writer.writerow({'first_name': 'Wonderful', 'last_name': 'Spam'})



def print_tally(counter):
    """
    tallies all of the strains in the world and outputs it.
    """
    ordered_bac_types = sorted(list(BACTERIA_TYPES))
    for bac_type in ordered_bac_types:
        line = "{}: {}".format(bac_type, counter[bac_type])
        print(line)

def count_bacteria_types(world):
    """
    takes a world and returns a counter of the bacteria types.
    """

    bac_strs = []
    for cell in world.cells:
        bac = cell.bacteria_type
        bac_str = bac.bacteria_type
        bac_strs.append(bac_str)
    counter = Counter(bac_strs)
    return counter

def map_world_to_csv(world, file_name):
    """
    creates map of world (location and bacterial cell type)
    """
    x_dimension = world.dimension_x_length
    y_dimension = world.dimension_y_length


    lines = []
    for y in range(y_dimension):
        line = []
        for x in range(x_dimension):
            cell = world.get_cell(x,y)
            bacteria_type_str = cell.bacteria_type.bacteria_type
            line.append(bacteria_type_str)
        lines.append(line)

    with open(file_name, "w") as file_handle:
        writer = csv.writer(file_handle)
        writer.writerows(lines)



def run_replicate(is_structured, number_of_generations, length_of_world,
    seed_proportions, file_name):
    """
    seed_proportions can take in any bacteria type and value between 0 and 1
    """
    earth = world.World(length_of_world, length_of_world,[])
    world.seed_world(earth, seed_proportions)


    counters = []
    counter = count_bacteria_types(earth)
    counters.append(counter)
    for gen_num in range(number_of_generations):
        earth.advance_one_generation(is_structured)
        counter = count_bacteria_types(earth)
        counters.append(counter)
        map_world_to_csv(earth, "Map_{}.csv".format(gen_num))

    output_counts_to_csv(file_name, counters)
    return earth

num_of_reps = 2
for rep_num in range(num_of_reps):
    filename = "teststructured_{}.csv".format(rep_num)
    World = run_replicate(is_structured=True, number_of_generations=1000,
        length_of_world=10, seed_proportions={"C":0.25, "S":0.25, "CL-": 0.25},
        file_name=filename)
