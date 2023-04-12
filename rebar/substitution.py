class Substitution:
    def __init__(self, substitution):
        self.substitution = substitution
        self.coord = int(substitution[1:-1])
        self.ref = substitution[0]
        self.alt = substitution[-1]

    def __repr__(self):
        return self.substitution

    def __lt__(self, other):
        return self.coord < other.coord

    def __le__(self, other):
        return self.coord <= other.coord

    def __gt__(self, other):
        return self.coord > other.coord

    def __ge__(self, other):
        return self.coord >= other.coord

    def __eq__(self, other):
        return self.substitution == other.substitution

    def __ne__(self, other):
        return self.substitution != other.substitution

    def __sub__(self, other):
        return self.coord + other.coord

    def __add__(self, other):
        return self.coord + other.coord

    def __hash__(self):
        return hash(str(self))
