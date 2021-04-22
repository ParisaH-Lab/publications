import os, shutil, sys
from translate_rotate import *


class pdbfile:
    '''
	read_file : readpdb file and returns list
	get_atoms    : takes a pdb list and specific atoms and returns them 
	replace_xyz  : takes a file and transform its coordinates
	'''

    def write_file(self, list_of_lines, filename='new.pdb'):
        fl = open(filename, 'w')
        for line in list_of_lines:
            fl.write(line)
        fl.close()


    # Requires floating point number
    # Returns floating point with correct
    # number of digits for pdb
    def set_number_digits(self, number):
        return '%.3f' % number


    # @Requires pdbfile
    # @Returns file as list
    def read_file(self, file):
        f = open(file, 'r')
        fl = f.readlines()
        f.close()
        return fl

    # @Requires pdblist
    # @Requires list with the desired atoms
    # @Requires residue - necessary to specify residue name and number
    def get_atoms(self, pdblist, list_of_atomnames, residuename):
        dict_of_coordinates = {}
        for line in pdblist:
            if str(line[17:20]).strip() == residuename:
                atomname = str(line[12:16]).strip()
                if atomname in list_of_atomnames:
                    coordinates = self.get_xyz(line)

                    dict_of_coordinates[atomname] = coordinates
        return dict_of_coordinates

    # @requires atom index and residue name
    def get_atoms_per_index(self, pdblist, list_of_atomindex, residuename):
        dict_of_coordinates = {}

        for line in pdblist:
            if str(line[17:20]).strip() == residuename:
                atomindex = str(line[7:11]).strip()

                if atomindex in list_of_atomindex:
                    coordinates = self.get_xyz(line)

                    dict_of_coordinates[atomindex] = coordinates

        return dict_of_coordinates


    # @Requires pdblist
    # @Requires list with the desired atoms
    # @Requires residue - necessary to specify residue name and number
    # which number is specified here
    def get_all_residueatoms(self, pdblist, residuename):

        dict_of_coordinates = {}

        for line in pdblist:

            if str(line[17:20]).strip() == residuename:
                # a fix needs to be inserted here
                if (1 == 2):
                    atomname = str(line[12:16]).strip()
                else:
                    atomname = str(line[7:11]).strip()

                coordinates = self.get_xyz(line)
                dict_of_coordinates[atomname] = coordinates

        return dict_of_coordinates

    # @Requires pdblist
    # @Requires list with the desired atoms
    # @Requires residue - necessary to specify residue name and number
    # which number is specified here
    def get_all_residueatoms(self, pdblist, residuename, residuenumber, dict_of_coordinates ):

        for line in pdblist:

            if str(line[17:20]).strip() == residuename and str(line[22:26]).strip() == residuenumber :
                if( line[13:15].strip() in ["CA","CB","N"]  ):
                    atomname = str(line[13:17]).strip()
                    coordinates = self.get_xyz(line)
                    # key for residue atoms:
                    tmp_key = residuename+"_"+atomname

                    dict_of_coordinates[ tmp_key ] = coordinates

        return dict_of_coordinates


    def get_protein_coordinates(self, pdblist):
        protein = []
        for line in pdblist:
            if line[0:4] == 'ATOM':
                protein.append(line)
        return protein


    # Requires a number
    # Set the right length for the number
    def set_length_digit(self, number):
        number = str(number)
        lngth = len(number)
        if lngth == 7:
            return ' ' + number
        if lngth == 6:
            return '  ' + number
        if lngth == 5:
            return '   ' + number
        if lngth == 4:
            return '    ' + number
        else:
            return number


    # Require file
    # Input crystal coordinates
    def get_transformed_coordinates(self, coordinates, pdbfile):
        t_coordinates = []
        for line in pdbfile:
            # fix needs to insert and if else here
            atomname = str(line[12:16]).strip()
            atomindex = str(line[7:11]).strip()
            atomname = atomindex
            if atomname in coordinates:
                x = '%.3f' % coordinates[atomname][0]
                y = '%.3f' % coordinates[atomname][1]
                z = '%.3f' % coordinates[atomname][2]

                x = self.set_length_digit(x)
                y = self.set_length_digit(y)
                z = self.set_length_digit(z)

                line = str(line[0:30]) + x + y + z + str(line[55:])
            t_coordinates.append(line)
        # Fix 18-03-2013
        t_coordinates.append("END\n")
        return t_coordinates


    def get_xyz(self, pdbline):
        x = float(pdbline[30:39])
        y = float(pdbline[39:47])
        z = float(pdbline[47:55])
        return array([x, y, z])


    def split_into_multiple_pdb_files(self, list_of_lines, delimiter='MODEL'):
        number = 0
        pdbfile = []
        for line in list_of_lines:

            if (line[0:5] == delimiter):
                if (number != 0):
                    filename = str(number) + '.pdb'
                    self.write_file(pdbfile, filename)
                number += 1
                pdbfile = []

            if (line[0:4] == 'ATOM'):
                pdbfile.append(line)



    # @Requires pdblist
    # @Requires list with the desired atoms
    # @Requires residue - necessary to specify residue name and number
    # which number is specified here
    def get_pdbcoordinates(self, pdblist ):
        dict_of_coordinates = {}
        for line in pdblist:

            if( line[0:4] == "ATOM" or line[0:4] == "HETA"):

                residuename = line[17:20].strip()
                residuenumber = str(line[22:26]).strip()
                atomname = str(line[13:17]).strip()

                coordinates = self.get_xyz(line)

                tmp_key = residuename+"_"+residuenumber+"_"+atomname

                dict_of_coordinates[ tmp_key ] = coordinates

        return dict_of_coordinates


    def get_transformed_coordinates_pdbfile(self, coordinates, pdbfile):
        t_coordinates = []
        for line in pdbfile:

            residuename = line[17:20].strip()
            residuenumber = str(line[22:26]).strip()
            atomname = str(line[12:16]).strip()

            tmp_key = residuename+'_'+residuenumber+'_'+atomname


            if tmp_key in coordinates:
                x = '%.3f' % coordinates[tmp_key][0]
                y = '%.3f' % coordinates[tmp_key][1]
                z = '%.3f' % coordinates[tmp_key][2]

                x = self.set_length_digit(x)
                y = self.set_length_digit(y)
                z = self.set_length_digit(z)

                line = str(line[0:30]) + x + y + z + str(line[55:])
                t_coordinates.append(line)
                # Fix 18-03-2013
                #t_coordinates.append("END\n")

        assert( len(t_coordinates) > 0)

        return t_coordinates


