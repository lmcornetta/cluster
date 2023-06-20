import molecule as mol
import pandas as pd
import sys, time

def molecules(filename):
    # Initial time
    start = time.time()

    # Transform the input text file into a pandas DataFrame
    df = pd.read_table(filename, sep='\s+')

    # Create a list of 'Molecule' instances read from given file (frame)
    list_of_molecules = []

    # Count indexes for totals
    number_of_ethanol, number_of_water = 0, 0

    # Begin iteration through the rows of the data frame
    # The loop reads each line (e.g. atom) and allocates it
    # in its original molecule.
    molecule = mol.Molecule()
    for index, row in df.iterrows():
        # Throw away dummy atoms from TIP4 force field
        if row['Atom'] == "MW":
            continue  
        # Read atom
        if molecule.is_complete() == False:
            molecule.include_atom(row['Atom'], float(row['x']), float(row['y']), float(row['z']))
        else:
            if molecule.name == "ethanol":
                number_of_ethanol += 1
            elif molecule.name == "water":
                number_of_water += 1
            molecule.center_of_mass()
            list_of_molecules.append(molecule)
            #molecule.print_molecule()
            del molecule
            molecule = mol.Molecule()
            molecule.include_atom(row['Atom'], float(row['x']), float(row['y']), float(row['z']))

    # Last molecule
    # Throw away dummy atoms from TIP4 force field
    if row['Atom'] == "MW":
        pass  
    # Read atom
    if molecule.is_complete() == False:
        molecule.include_atom(row['Atom'], float(row['x']), float(row['y']), float(row['z']))
    else:
        if molecule.name == "ethanol":
            number_of_ethanol += 1
        elif molecule.name == "water":
            number_of_water += 1
        molecule.center_of_mass()
        list_of_molecules.append(molecule)
        #molecule.print_molecule()
        del molecule

    end = time.time()

    # List done!
    print(" ==> List of Molecules built! <== ")
    print(" -- elapsed time to generate list of objects: %f sec.\n"%(end - start))

    # Print status of found molecules
    print(" ==> Totals <== ")
    print("Total number of water molecules: " + str(number_of_water))
    print("Total number of ethanol molecules: " + str(number_of_ethanol))
    print(" -- ratio #ethanol/#total (ethanol concentration): %f"%((number_of_ethanol)/(number_of_ethanol + number_of_water)))
    return list_of_molecules, number_of_ethanol, number_of_water
