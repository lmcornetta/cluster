import sys
from read_frame import molecules

if __name__=="__main__":
    # Lets assume we have the list of molecules generated by the read_frame.py
    filename = sys.argv[1]
    list_molecules, number_of_ethanol, number_of_water = molecules(filename)
    total_molecules = number_of_ethanol + number_of_water

    # Choose number of clusters to generate
    number_of_clusters = 120

    # Choose radius cutoff, in angstroms
    radius_cutoff = 8.40

    # Choose central molecule {options: "ethanol" or "water"}
    central_molecule_type = "ethanol"

    cluster, index = 1, 0
    surf, bulk = 1, 1
    surf_ymin, surf_ymax = 30.0, 50.0
    while cluster <= number_of_clusters:

        # Central molecule
        central_molecule = list_molecules[index]
        if central_molecule.name != central_molecule_type or central_molecule.is_broken():
            index += 1
            continue

        # Output file
        if central_molecule.ycm > surf_ymin and central_molecule.ycm < surf_ymax:
            fcluster = open("surf%d.xyz"%surf,"w")
            fmolcas = open("surf%d.input"%surf,"w")
            surf += 1
        else:
            fcluster = open("bulk%d.xyz"%bulk,"w")
            fmolcas = open("bulk%d.input"%bulk,"w")
            bulk += 1

        fmolcas.write("&SEWARD\n Cholesky\n")
        h_index = 1
        c_index = 1
        o_index = 1
        for i in range(central_molecule.number_of_atoms):
            if central_molecule.list_of_types[i] == 'H':
                fmolcas.write(" Basis set\n  H.ANO-RCC-VDZP\n")
                fmolcas.write(" H" + str(h_index) + "\t" + str(central_molecule.list_of_coord[i][0]) + "\t" + str(central_molecule.list_of_coord[i][1]) +
                              "\t" + str(central_molecule.list_of_coord[i][2]) + " Angstrom" + "\n")
                fmolcas.write(" End of basis\n")
                h_index += 1

            if central_molecule.list_of_types[i] == 'C':
                fmolcas.write(" Basis set\n  C.ANO-RCC-VDZP\n")
                fmolcas.write(" C" + str(c_index) + "\t" + str(central_molecule.list_of_coord[i][0]) + "\t" + str(central_molecule.list_of_coord[i][1]) +
                              "\t" + str(central_molecule.list_of_coord[i][2]) + " Angstrom" + "\n")
                fmolcas.write(" End of basis\n")
                c_index += 1

            if central_molecule.list_of_types[i] == 'O':
                fmolcas.write(" Basis set\n  O.ANO-RCC-VDZP\n")
                fmolcas.write(" O" + str(o_index) + "\t" + str(central_molecule.list_of_coord[i][0]) + "\t" + str(central_molecule.list_of_coord[i][1]) +
                              "\t" + str(central_molecule.list_of_coord[i][2]) + " Angstrom" + "\n")
                fmolcas.write(" End of basis\n")
                o_index += 1

        # All good! Let's go
        molecules_in_cluster = [central_molecule]
        cluster_size = 1
        cluster_atoms_quantity = central_molecule.number_of_atoms
        for molecule in list_molecules:
            if molecule.is_broken():
                continue
            d = central_molecule.distance_to_other(molecule)
            if d < radius_cutoff and d > 0.0:
                molecules_in_cluster.append(molecule)
                # Writting to molcas input
                for i in range(molecule.number_of_atoms):
                    if molecule.list_of_types[i] == 'H':
                        fmolcas.write(" Basis set\n  H.ANO-S-MB\n")
                        fmolcas.write(" H" + str(h_index) + "\t" + str(molecule.list_of_coord[i][0]) + "\t" + str(molecule.list_of_coord[i][1]) +
                                      "\t" + str(molecule.list_of_coord[i][2]) + " Angstrom" + "\n")
                        fmolcas.write(" End of basis\n")
                        h_index += 1

                    if molecule.list_of_types[i] == 'C':
                        fmolcas.write(" Basis set\n  C.ANO-S-MB\n")
                        fmolcas.write(" C" + str(c_index) + "\t" + str(molecule.list_of_coord[i][0]) + "\t" + str(molecule.list_of_coord[i][1]) +
                                      "\t" + str(molecule.list_of_coord[i][2]) + " Angstrom" + "\n")
                        fmolcas.write(" End of basis\n")
                        c_index += 1

                    if molecule.list_of_types[i] == 'O':
                        fmolcas.write(" Basis set\n  O.ANO-S-MB\n")
                        fmolcas.write(" O" + str(o_index) + "\t" + str(molecule.list_of_coord[i][0]) + "\t" + str(molecule.list_of_coord[i][1]) +
                                      "\t" + str(molecule.list_of_coord[i][2]) + " Angstrom" + "\n")
                        fmolcas.write(" End of basis\n")
                        o_index += 1
            
                cluster_size += 1
                cluster_atoms_quantity += molecule.number_of_atoms
        fmolcas.write(" Douglas-Kroll\n&SCF\n >>> COPY $Project.ScfOrb $CurrDir/$Project.ScfOrb")
        fmolcas.close()
        fcluster.write("%d\n\n"%cluster_atoms_quantity)
        for i in range(cluster_size):
            molecules_in_cluster[i].write_molecule(fcluster)
        fcluster.close()
        cluster += 1
        index += 1
    
    # Print totals
    print("\n Total number of generated clusters with radius cutoff of %f angstroms: %d"%(radius_cutoff, cluster - 1))
    print(" ---> Bulk clusters: %d"%(bulk - 1))
    print(" ---> Surface clusters: %d\n"%(surf - 1))
