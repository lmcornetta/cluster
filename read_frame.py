import molecule as mol
import pandas as pd
import matplotlib.pyplot as plt
import time
import numpy as np

def molecules(filename):

    # Global parameters for the plots
    plt.rcParams['font.family'] = 'arial'
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['legend.fontsize'] = 17
    plt.rcParams['xtick.labelsize'] = 17
    plt.rcParams['ytick.labelsize'] = 17
    
    # Colormap. These colors are recomended for publications where a 'greyscale' 
    # version is required 
    new_blue = '#004488'
    new_yellow = '#DDAA33'
    new_red = '#BB5566'

    # Initial time
    start = time.time()

    # Transform the input text file into a pandas DataFrame
    df = pd.read_table(filename, sep='\s+')

    # Plot histograms of atomic coordinates. Just for checking.
    # Sometimes the slab is located at the edge of the box, sometimes it is
    # in the middle... Plotting the histograms is a good practive to understand
    # the arrangement. Uncomment below if needed.
    ### plt.hist(df['x'])
    ### plt.show()
    ### plt.hist(df['y'])
    ### plt.show()
    ### plt.hist(df['z'])
    ### plt.show()

    # Create a list of 'Molecule' instances read from given file (frame)
    list_of_molecules = []
    list_of_ethanol_ycm = []
    list_of_water_ycm = []

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
            molecule.center_of_mass()
            if molecule.name == "ethanol":
                number_of_ethanol += 1
                list_of_ethanol_ycm.append(molecule.ycm)
            elif molecule.name == "water":
                number_of_water += 1
                list_of_water_ycm.append(molecule.ycm)
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
        molecule.center_of_mass()
        if molecule.name == "ethanol":
            number_of_ethanol += 1
            list_of_ethanol_ycm.append(molecule.ycm)
        elif molecule.name == "water":
            number_of_water += 1
            list_of_water_ycm.append(molecule.ycm)
        list_of_molecules.append(molecule)
        #molecule.print_molecule()
        del molecule

    end = time.time()

    # List done!
    print(" ==> List of Molecules built! <== ")
    print(" -- elapsed time to generate list of objects: %f sec.\n"%(end - start))
    ratio = (number_of_ethanol)/(number_of_ethanol + number_of_water)

    # Distribution and densities (figure)
    fig, ax = plt.subplots()

    # Calculating densities
    count_ethanol_, bins_ethanol_ = np.histogram(list_of_ethanol_ycm, bins=200)
    count_water_, bins_water_ = np.histogram(list_of_water_ycm, bins=bins_ethanol_)
    ax.plot(bins_ethanol_[:-1], 7*count_ethanol_, ls='-', color=new_yellow, lw=2.5)
    ax.plot(bins_water_[:-1], 7*count_water_, ls='-', color=new_red, lw=2.5)
    ax.set_xlim(0.0, 150.0)

    ax_ = ax.twinx()
    ax_.plot(bins_ethanol_[:-151], np.array(count_ethanol_[:-150])/(np.array(count_ethanol_[:-150]) + np.array(count_water_[:-150])), ls='-', lw=2.8, color=new_blue)
    ax_.plot(bins_ethanol_[:-151], ratio*np.ones(50), color=new_blue, lw=2.0, ls='--')
    fig.text(0.41, 1.6*ratio, '{:.2f}'.format(ratio), fontsize=19, color=new_blue)
    ax_.set_ylabel("ethanol concentration", fontsize=19, rotation=270, labelpad=20, color=new_blue)
    ax_.set_ylim(0.3*ratio, 1.0)
    ax_.tick_params(axis='y', labelcolor=new_blue)

    # Plot histograms of CMs along the y-direction
    counts_ethanol, bins_ethanol = np.histogram(list_of_ethanol_ycm)
    ax.stairs(counts_ethanol, bins_ethanol, color=new_yellow, lw=1.8)

    counts_water, bins_water = np.histogram(list_of_water_ycm)
    ax.stairs(counts_water, bins_water, color=new_red, lw=1.8)
    
    ax.hist(list_of_ethanol_ycm, alpha=0.5, color=new_yellow, label='ethanol')
    ax.hist(list_of_water_ycm, alpha=0.5, color=new_red, label='water')
    ax.set_title(r"Distribution of CMs along the $z-$axis", fontsize=18)
    ax.set_xlabel(r"$z$ coordinate (angs.)", fontsize=18)
    ax.set_ylabel("counts", fontsize=18)
    ax.legend(loc='upper right')
    
    fig.tight_layout()
    plt.show()

    # Print status of found molecules
    print(" ==> Totals <== ")
    print("Total number of water molecules: " + str(number_of_water))
    print("Total number of ethanol molecules: " + str(number_of_ethanol))
    print(" -- ratio #ethanol/#total (ethanol concentration): %f"%ratio)
    return list_of_molecules, number_of_ethanol, number_of_water
