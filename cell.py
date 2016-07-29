
BACTERIA_TYPES = {"C", "L+", "L-", "CL+", "CL-", "S", "E"}

COST_FOR_CARRYING_COLICIN = 0.05
COST_FOR_CARRYING_PHAGE = 0.03

class Bacteria:
    """
    A class that represent the bacterial genotypes in our model
    (C, L+, L-,  CL+, CL-, S, E)
    """
    def __init__(self, bacteria_type):
        assert bacteria_type in BACTERIA_TYPES
        self.bacteria_type = bacteria_type

    def __str__(self):
        return self.bacteria_type

    def is_immune_to_colicin(self):
        return "C" in self.bacteria_type

    def is_immune_to_phage(self):
        return "L" in self.bacteria_type

    def can_release_colicin(self):
        return self.bacteria_type in {"CL+", "CL-"}

    def can_release_phage(self):
        return "L+" in self.bacteria_type

    def can_lyse(self):
        return self.can_release_colicin() or self.can_release_phage() or self.bacteria_type == "L-"

    def can_become_lysogen(self):
        return self.bacteria_type in {"C", "S"}

    def convert_to_lysogen(self):
        assert self.can_become_lysogen()
        if self.bacteria_type == "S":
            return Bacteria("L+")
        if self.bacteria_type == "C":
            return Bacteria("CL+")

    def replication_rate(self):
        """
        returns the probability of sucessful replication.
        E = 0
        S = 1
        C = 1 - COST_FOR_CARRYING_COLICIN
        ...
        """
        if self.bacteria_type == "E":
            return 0
        if self.bacteria_type == "S":
            return 1
        rate = 1.0
        if "C" in self.bacteria_type:
            rate = rate - COST_FOR_CARRYING_COLICIN
        if "L" in self.bacteria_type:
            rate = rate - COST_FOR_CARRYING_PHAGE
        return rate


class Cell:
    """
    A cell represents a position in a grid. The cell can have various things in
    it, including colicin, and or phage, and or bacteria, or be empty.
    """
    def __init__(self, has_colicin=False, has_phage=False, bacteria_type=Bacteria("E")):
        self.has_colicin = has_colicin
        self.has_phage = has_phage
        self.bacteria_type = bacteria_type

    def __str__(self):
        return "Cell(has_colicin={}, has_phage={}, bacteria_type={})".format(
            self.has_colicin, self.has_phage,
            self.bacteria_type)

    def __repr__(self):
        return str(self)
