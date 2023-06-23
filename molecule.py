import numpy as np

# Defining atomic labels in the simulation.
# This is important for identifiyng which molecule we are dealing with
ethanol_labels = ["CA", "CB", "HA1", "HA2", "HB1", "HB2", "HB3", "OH", "HO"]
water_labels = ["OW", "HW1", "HW2"]
dummy_labels = ["MW"]

carbon_labels = ["CA", "CB"]
oxygen_labels = ["OH", "OW"]
hydrogen_labels = ["HA1", "HA2", "HB1", "HB2", "HB3", "HW1", "HW2", "HO"]

# Dictionary for atomic masses
atomic_mass = {'CA' : 12.011,
               'CB' : 12.011,
               'OH' : 15.999,
               'OW' : 15.999,
               'HA1': 1.0078,
               'HA2': 1.0078,
               'HB1': 1.0078,
               'HB2': 1.0078,
               'HB3': 1.0078,
               'HW1': 1.0078,
               'HW2': 1.0078,
               'HO' : 1.0078}

# Class Molecule
class Molecule:

    def __init__(self):
        # Initialization (empty Molecule)
        self.number_of_atoms = 0
        self.name = None
        self.list_of_atoms = []
        self.list_of_types = []
        self.list_of_coord = []

    def is_complete(self):
        # Check if there is no further atoms to insert
        # in order for the Molecule to be complete
        if self.name == "ethanol" and self.number_of_atoms == 9:
            return True
        elif self.name == "water" and self.number_of_atoms == 3:
            return True
        else:
            return False
    
    def include_atom(self, label, x, y, z):
        # Check if this is the first atom to be
        # inserted in the Molecule
        if self.name == None:
            if label in ethanol_labels:
                self.name = 'ethanol'
            if label in water_labels:
                self.name = 'water'

        self.list_of_atoms.append(label)
        if label in carbon_labels:
            self.list_of_types.append("C")
        if label in oxygen_labels:
            self.list_of_types.append("O")
        if label in hydrogen_labels:
            self.list_of_types.append("H")
        self.list_of_coord.append(np.array([x,y,z], dtype=np.float32))
        self.number_of_atoms += 1

    def center_of_mass(self):
        # Check it there are more atoms to fill
        if self.is_complete() == False:
            raise Exception("Molecule is not complete. Can't calculate CM.")

        # Calculate the CM coordinates of the Molecule
        xcm = 0.0
        ycm = 0.0
        zcm = 0.0
        mol_mass = 0.0
        for i in range(self.number_of_atoms):
            mi = atomic_mass[self.list_of_atoms[i]]
            mol_mass += mi
            xcm += mi*self.list_of_coord[i][0]
            ycm += mi*self.list_of_coord[i][1]
            zcm += mi*self.list_of_coord[i][2]
        self.xcm = xcm/mol_mass
        self.ycm = ycm/mol_mass
        self.zcm = zcm/mol_mass
    
    def is_broken(self):
        # Check if the molecule is broken. That means it is disposed
        # partially at one part of the box and partially at another.
        # That happens sometimes due to periodic conditions.
        for i in range(self.number_of_atoms):
            dx = np.abs(self.xcm - self.list_of_coord[i][0])
            dy = np.abs(self.ycm - self.list_of_coord[i][1])
            dz = np.abs(self.zcm - self.list_of_coord[i][2])
            if dx > 5.0 or dy > 5.0 or dz > 5.0:
                return True
        return False

    def distance_to_other(self, other):
        # Calculate distance to other Molecule CM, in angstroms
        # This molecule
        xcm1 = self.xcm
        ycm1 = self.ycm
        zcm1 = self.zcm

        # The "other" molecule
        xcm2 = other.xcm
        ycm2 = other.ycm
        zcm2 = other.zcm

        # Cartesian distance
        d12 = np.sqrt((xcm1 - xcm2)**2 + (ycm1 - ycm2)**2 + (zcm1 - zcm2)**2)
        return d12
    
    def write_molecule(self, fbuffer):
        # Write molecule to xyz file
        for i in range(self.number_of_atoms):
            fbuffer.write(str(self.list_of_types[i]) + "\t" + str(self.list_of_coord[i][0]) + "\t" + str(self.list_of_coord[i][1]) + "\t" + str(self.list_of_coord[i][2]) + "\n")

    def print_molecule(self):
        print(" ==> Molecule found <==")
        print("   Name: %s"%self.name)
        print("   Number of atoms: %s"%self.number_of_atoms)
        print(self.list_of_atoms)
        print(self.list_of_types)
        print(self.list_of_coord)
        print("Center of mass: " + str(self.center_of_mass()) + "\n")
