from cell import Bacteria, Cell
import random


PROB_OF_INDUCTION = 0.0037
COLICIN_PERSISTENCE =  10
PHAGE_PERSISTENCE = 3

def seed_world(world, seed_proportions):
    total_proportion = sum(seed_proportions.values())
    assert total_proportion <= 1, "Seed proportions must tally to less than or equal 1"


    number_of_positions = len(world.cells)
    new_cells = []
    for bacteria_type_str, proportion in seed_proportions.items():
        num_cells_of_type = int(number_of_positions * proportion)
        for _ in range(num_cells_of_type):
            bacteria_type = Bacteria(bacteria_type_str)
            new_cells.append(Cell(bacteria_type=bacteria_type))

    while len(new_cells) < number_of_positions:
        new_cells.append(Cell())
    random.shuffle(new_cells)
    world.cells = new_cells

class World:
    """
    A grid of cells.
    """
    def __init__(self, dimension_x_length, dimension_y_length, cells=[]):
        self.dimension_x_length = dimension_x_length
        self.dimension_y_length = dimension_y_length
        assert self.dimension_x_length >= 1
        assert self.dimension_y_length >= 1
        self.cells = cells[:]
        total_cell_count = self.dimension_x_length * self.dimension_y_length
        assert len(self.cells) <= total_cell_count
        while len(self.cells) < total_cell_count:
            self.cells.append(Cell())

    def get_all_coordinates_in_world(self):
        """
        Go through all coordinates by row, starting with (0,0), (1,0), (2,0)
        """
        results = []
        for y in range(self.dimension_y_length):
            for x in range(self.dimension_x_length):
                coord = (x,y)
                results.append(coord)
        return results


    def advance_one_generation(self, is_structured):
        """
        Colicin releasers relase colicin.
        Lysogens with functional phage release phage.
        Cells without immunity to colicin or phage will die if either is present.
        Remove free colicin and free phage.
        All alive bacteria replicate once.
        """
        self.prob_induce_lysis_of_all_cells(PROB_OF_INDUCTION)
        self.kill_non_immunes()
        self.degrade_colicin_and_phage()
        self.replicate(is_structured)



    def induce_lysis(self, x, y):
        """
        Induces lysis and releases colicin and/or phage into nearest neighbors
        set by radius for bacterial cells that can lyse ( anything with L )
        """
        cell = self.get_cell(x,y)
        bac = cell.bacteria_type
        release_phage = bac.can_release_phage()
        release_colicin = bac.can_release_colicin()

        nearest_neighbors = self.get_nearest_neighbors_of_cell(x, y, radius=2)
        for neighbor in nearest_neighbors:
            if release_phage:
                cell.has_phage = PHAGE_PERSISTENCE
                neighbor.has_phage = PHAGE_PERSISTENCE
            if release_colicin:
                cell.has_colicin = COLICIN_PERSISTENCE
                neighbor.has_colicin = COLICIN_PERSISTENCE
        if bac.can_lyse():
            cell.bacteria_type = Bacteria("E")




    def prob_induce_lysis_of_all_cells(self, induction_rate):
        """
        for each cell in the world, the cell will be induced some proportion
        of the time randomly. If an induced cell can produce colicin or phage,
        such will be added to its neighbors.
        """
        all_coords = self.get_all_coordinates_in_world()
        for coord in all_coords:
            x,y = coord
            induced = random.random() < induction_rate
            if induced:
                self.induce_lysis(x, y)



    def kill_non_immunes(self):
        """
        kills bacterial cells that are sensitive to colicin or phage if colicin
        or phage is present in cell.
        """
        all_coords = self.get_all_coordinates_in_world()
        for coord in all_coords:
            x,y = coord
            cell = self.get_cell(x,y)
            bac = cell.bacteria_type
            is_immune_to_colicin = bac.is_immune_to_colicin()
            if not is_immune_to_colicin and cell.has_colicin:
                cell.bacteria_type = Bacteria("E")
            is_immune_to_phage = bac.is_immune_to_phage()
            if not is_immune_to_phage and cell.has_phage:
                cell.bacteria_type = Bacteria("E")

    def degrade_colicin_and_phage(self):
        """
        removes all colicin and phage from cells
        """
        all_coords = self.get_all_coordinates_in_world()
        for coord in all_coords:
            x,y = coord
            cell = self.get_cell(x,y)
            if cell.has_colicin:
                cell.has_colicin -= 1
            if cell.has_phage:
                cell.has_phage -= 1

    def replicate(self, is_structured):
        """
        All non E bacteria will replicate randomly into a neighbor.
        """
        all_coords = self.get_all_coordinates_in_world()
        random.shuffle(all_coords)
        for coord in all_coords:
            x,y = coord
            cell = self.get_cell(x,y)
            bac = cell.bacteria_type
            replication_prob = bac.replication_rate()
            if random.random() < replication_prob:
                if is_structured:
                    four_neighbors = self.get_nearest_four_neighbors_of_cell(x,y)
                else:
                    four_neighbors = self.get_random_neighbors_of_cell(x, y, 4)
                unlucky = random.choice(four_neighbors)
                unlucky.bacteria_type = Bacteria(bac.bacteria_type)


    def check_bounds(self, x, y):
        """
        checks the bounds of the world
        """
        assert x >= 0
        assert x < self.dimension_x_length
        assert y >= 0
        assert y < self.dimension_y_length

    def __str__(self):
        """
        returns a string
        """
        return "World(dimension_x_length={}, dimension_y_length={}, cells={})".format(
        self.dimension_x_length, self.dimension_y_length, self.cells)

    def get_cell(self, x, y):
        """
        this will return a cell associated with the coordinates x,y. indexed by 0.
        """
        self.check_bounds(x,y)
        index = x + y * self.dimension_x_length
        return self.cells[index]

    def get_nearest_four_neighbors_of_cell(self, x, y):
        """
        gets nearest four neighbors of a focal cell: in a grid this is the cell
        immediately above, below, to the right, and left of a focal cell.
        """
        self.check_bounds(x,y)
        result = []
        result.append(self.get_wrapped_cell(x, y-1))
        result.append(self.get_wrapped_cell(x+1, y))
        result.append(self.get_wrapped_cell(x, y+1))
        result.append(self.get_wrapped_cell(x-1, y))
        return result

    def get_nearest_neighbors_of_cell(self, x, y, radius):
        """
        Gets all the cells within a radius of the focal cell.
        Excludes the focal cells.

        Radius of 0 is nothing.
        Radius of 1 is the nearest 8 neighbors.
        Radius of 2 is the nearest 24 neighbors.
        """
        self.check_bounds(x,y)
        x_range = list(range(x - radius, x + radius + 1))
        y_range = list(range(y - radius, y + radius + 1))
        result = []
        for x_coord in x_range:
            for y_coord in y_range:
                if x_coord == x and y_coord == y:
                    continue
                result.append(self.get_wrapped_cell(x_coord, y_coord))
        return result



    def get_random_neighbors_of_cell(self, x, y, number_of_neighbors):
        """
        gets random cells from world not including focal cell.
        for well-mixed environment.
        """
        self.check_bounds(x,y)
        focal_coordinate = (x,y)
        coordinates = set()
        while len(coordinates) < number_of_neighbors:
            rand_x = random.randrange(self.dimension_x_length)
            rand_y = random.randrange(self.dimension_y_length)
            rand_coord = (rand_x, rand_y)
            if rand_coord == focal_coordinate:
                continue
            coordinates.add(rand_coord)
        results = []
        for coord in coordinates:
            x, y = coord
            rand_cell = self.get_cell(x, y)
            results.append(rand_cell)
        return results




    def get_wrapped_cell(self, x, y):
        """
        If given a coordinate outside the bounds, this will return the cell after
        wrapping around the boundaries.
        """
        wrapped_x = x % self.dimension_x_length
        wrapped_y = y % self.dimension_y_length
        return self.get_cell(wrapped_x, wrapped_y)
