#!/usr/bin/python

import numpy as np
import fractions
import math as m
import os
import decimal as dc
import sys
import json
import time

#
# Material Builder 
# Author: Jakub Kaminski, UCLA 2015
#	Modified by Boya Song, UCLA Winter 2015
#
# Based on Interface Builder
#
if len(sys.argv) <=1:
	print "Usage:\n%s optionsFile"%sys.argv[0]
	exit()

inputFile = sys.argv[1]

def createMillerList(maxMillerInd):
	 # create a list that contains strings for all the possible Miller indices
	 # from "0 0 1" to "maxMillerInd maxMillerInd maxMillerInd"

	MillerList=list() 
	x1=0
	x2=0
	x3=0

	while x1 <= maxMillerInd:
		x2=0
		while x2 <= maxMillerInd:
			x3=0
			while x3 <= maxMillerInd:
				MillerList.append(str(x1)+" "+str(x2)+" "+str(x3))
				x3 = x3+1
			x2=x2+1
		x1=x1+1

	MillerList.remove("0 0 0") # we cannot have "0 0 0"

	return MillerList

def readInput(inputFile):

	file = open(inputFile,'r')
	line = file.readline()
	line = line.split()
	subCIF = line[1]

	line = file.readline()
	line = line.split()
	useMillerList = bool(int(line[1]))

	line = file.readline()
	line = line.split()
	maxMillerInd = int(line[1])

	line = file.readline()
	line = line.split()
	MillerIndList = " ".join(line[1:-7])
	MillerIndList = MillerIndList.split('|')

	line = file.readline()
	line = line.split()
	nL = float(line[1])

	line = file.readline()
    line = line.split()
    checkPolarity = bool(line[1])

    line = file.readline()
    line = line.split()
    z_digit = int(line[1])

    line = file.readline()
    line = line.split()
    dis_tol_rate = float(line[1])

    return subCIF, useMillerList, maxMillerInd, MillerIndList,\
     nL, checkPolarity, z_digit, dis_tol_rate



def getMillerFromString(millerString):

	out = np.array((0,0,0))
	millerString = millerString.split();
	for item in range(len(millerString)):
		out[item] = millerString[item]
	
	return out

def uqLabels(labels,types):
        # Find unique labels in the list
        # Return dictionary with key=integer number, value=labels
        tmp = []
	# Get all the unique atoms
        for i in labels:
                if i not in tmp:
                        tmp.append(i)

	# Get all the atoms that are in dictionary to the list
	vals = types.values()

	nextIdx = 0
	if types != {}:
		maxTyp = max(types)
		nextIdx = maxTyp + 1

	# Put the atoms to the dictionary
	for i in tmp:
		if i not in vals:
			types[nextIdx] = i
			nextIdx += 1

        return types

def ReadCIF(filename,atomTypes):
	# Reads information from CIF file:
	# - ICSD database code for the material
	# - unit cell lengths and angles
	# - evaluates xyz coordinates based on fractional coordinates and 
	# symmetry position of unique atoms
	# Returns:
	# - ICSD database code for the material
	# - unit cell lengths and angles
	# - list with atom labels
	# - array with xyz coordinates

	file = open(filename,'r')

	lines = file.readlines()
	file.close()

	# Strip whitespaces and line endings
	for i in range(len(lines)):
		lines[i]=lines[i].strip()

	for line in lines:
		if line.find("data") >= 0:
			MatID = line

		elif line.find("_cell_length_a") >= 0:
			line = line.split()
			cellA =  CheckParentheses(line[1])

		elif line.find("_cell_length_b") >= 0:
			line = line.split()
			cellB =  CheckParentheses(line[1])

		elif line.find("_cell_length_c") >= 0:
			line = line.split()
			cellC =  CheckParentheses(line[1])

		elif line.find("_cell_angle_alpha") >= 0:
			line = line.split()
			alpha = CheckParentheses(line[1])

		elif line.find("_cell_angle_beta") >= 0:
			line = line.split()
			beta = CheckParentheses(line[1])

		elif line.find("_cell_angle_gamma") >= 0:
			line = line.split()
			gamma = CheckParentheses(line[1])
		elif line.find("_symmetry_space_group_name_H-M") >= 0:
			SGsymbol = extractSpaceGroup(line)

	f2cmat = FracToCart(cellA,cellB,cellC,alpha,beta,gamma)
	print "SPACE GROUP: ",SGsymbol

	# READ ATOM SYMMETRY POSITION CORRESPONDING TO SPACE GROUP
	sgInput = open('spacegroups.json','r')
	sgData = json.load(sgInput)
	try:
		atompos = sgData[SGsymbol]
	except KeyError:
		print "Unrecognized spacegroup symbol read from CIF file %s"%SGsymbol
		exit()

	# READ ATOM FRACTIONAL COORDINATES
	# Find index of the the beginning of atom fractional coordinates
	# Find either _atom_site_label or _atom_site_symbol and see what is first
	# Set the ib to the fist one
	try:
		ibL = lines.index("_atom_site_label")
	except ValueError:
		ibL = 100000 # Just put big unresonable value

	try:
		ibS = lines.index("_atom_site_type_symbol")
	except ValueError:
		ibS = 100000 # Just put big unresonable value

	# Find where is the begining of atom definition loop in CIF file.
	if ibL < ibS:
		ib = ibL
	else:
		ib = ibS

	allLabelsRead = False
	ibSave = ib
	while not allLabelsRead:
		ib += 1
		if lines[ib][0] == "_":
			if lines[ib][0:18] == "_atom_site_fract_x":
				xPos = ib - ibSave
			if lines[ib][0:18] == "_atom_site_fract_y":
				yPos = ib - ibSave
			if lines[ib][0:18] == "_atom_site_fract_z":
				zPos = ib - ibSave
		else:
			allLabelsRead = True

	# Get atom positions
	# We don't know how many atoms ther is the structure
	# Read them unitl keyworld "loop", or line begining with "_"
	# or line begining with "#" or end of file is encountered

	frapos = []
	while 1:
		try:
			line = lines[ib]
		except IndexError:
			break # end of file reached 

		# If file finishes with empty line
		try:
			test = line[0]
		except IndexError:
			break

		if (line[0:4] != "loop") and (line[0] != "_") \
		   and (line[0] != "#"):
			frapos.append(line)
			ib += 1
		else:
			break

	# Now, knowing fractional coordinates and the symmetry positions, 
	# build fractional xyz coordinates.
	# First we will put fractional coordinates into array fracoord

	fracoord = np.zeros((len(frapos),3))
	atomlabels = []
	ii = 0 #index
	for atom in frapos:
		atom = atom.split()
		symbol = atom[0]
		# Strip the digits from the symbol name:
		if symbol[-1].isdigit():
			symbol = symbol[0:-1]

		x = float(CheckParentheses(atom[xPos]))
		y = float(CheckParentheses(atom[yPos]))
		z = float(CheckParentheses(atom[zPos]))

		fracoord[ii,0] = x
		fracoord[ii,1] = y
		fracoord[ii,2] = z
		atomlabels.append(symbol)
		ii = ii+1
	# Create array with fractional coordinates based on
	# symmetry positions for each atom type
#	positions = np.zeros((len(fracoord)*len(atompos),3))

	iatlab = 0 # index of symmetry unique atom
	atoms = []# list to be fileld with atom indexes after symmetry operation
	for atom in fracoord:
		postmp = []
		for pos in atompos: #For each element in atompos, convert it 
			            # to fraction, multiply by fractional 
				    # coordinate and save to positions array
#			pos = pos.split()
			pos = pos.split(",") # For group
			coordtmp = []
#			for elem in pos[1:]:
			for elem in pos:     # For group
				# Read symmetry operation and constuct 
				# fractional coordiante

				s,o = ReadSym(elem)

				coord = s[0]*o[0]*atom[0] + \
					s[1]*o[1]*atom[1] + \
					s[2]*o[2]*atom[2] + \
					s[3]*o[3]             

				# Normalise to [0,1)
				if coord < 0 : coord = 1+coord 
				if coord > 1 : coord = coord-1
# WARNING
				# Round to 4 digits to avoid rounding error
# END WARNING
				coord = round(coord,3)
				coordtmp.append(coord)
			#Find duplicates	
			if coordtmp not in postmp:
				postmp.append(coordtmp)
				atoms.append(atomlabels[iatlab])

		# If first atom or only one atom, create position array
		if iatlab == 0:
			positions = np.array(postmp)				
		else: # add 2nd and other atoms to position array
			tmp = np.array(postmp)
			positions = np.vstack([positions,tmp])
		
	# Find and add 1 to every 0 in postmp. Add new coordiante to positions
		for elem in postmp:
			if elem[0] == 0:
				newrow = np.array(elem)
				newrow[0] += 1
				diff = CmpRows(positions,newrow)
				if diff:
					positions = \
						np.vstack([positions,newrow])
					atoms.append(atomlabels[iatlab])
			if elem[1] == 0:
				newrow = np.array(elem)
				newrow[1] += 1
				diff = CmpRows(positions,newrow)
				if diff:
					positions = \
						np.vstack([positions,newrow])
					atoms.append(atomlabels[iatlab])
			if elem[2] == 0:
				newrow = np.array(elem)
				newrow[2] += 1
				diff = CmpRows(positions,newrow)
				if diff:
					positions = \
						np.vstack([positions,newrow])
					atoms.append(atomlabels[iatlab])
			if elem[0] == 0 and elem[1] == 0:
				newrow = np.array(elem)
				newrow[0] += 1
				newrow[1] += 1
				diff = CmpRows(positions,newrow)
				if diff:
					positions = \
						np.vstack([positions,newrow])
					atoms.append(atomlabels[iatlab])
			if elem[0] == 0 and elem[2] == 0:
				newrow = np.array(elem)
				newrow[0] += 1
				newrow[2] += 1
				diff = CmpRows(positions,newrow)
				if diff:
					positions = \
						np.vstack([positions,newrow])
					atoms.append(atomlabels[iatlab])
			if elem[1] == 0 and elem[2] == 0:
				newrow = np.array(elem)
				newrow[1] += 1
				newrow[2] += 1
				diff = CmpRows(positions,newrow)
				if diff:
					positions = \
						np.vstack([positions,newrow])
					atoms.append(atomlabels[iatlab])
			if elem[1] == 0 and elem[1] == 0 and elem[2] == 0:
				newrow = np.array(elem)
				newrow[0] += 1
				newrow[1] += 1
				newrow[2] += 1
				diff = CmpRows(positions,newrow)
				if diff:
					positions = \
						np.vstack([positions,newrow])
				atoms.append(atomlabels[iatlab])
		iatlab += 1 # next atom

        #Convert atoms list to numpy array with atom types. Atom types will be hold in dictonary
        # Get the dictionary
        atomTypes = uqLabels(atoms,atomTypes)
        # Convert atoms list to array
        ii = 0
        atomstmp = np.zeros(len(atoms))
        for label in atoms:
                for uqLabIdx,uqLab in atomTypes.iteritems():
                        if label == uqLab:
                                atomstmp[ii] = uqLabIdx
                ii+=1
        atoms = atomstmp

	positions = np.dot(positions,f2cmat)
	#TODO: check if this is always true.
	# This is done for the cases when if CIF file the coordinates of first atom are different from 0.0 0.0 0.0
	# Such cases prevents proper rotations on def plane, where the desired plane doesnt have z=0
	shift = positions[0].copy()
	positions -= shift
	#end TODO
	return MatID,f2cmat,atoms,positions,atomTypes
def CmpRows(mat1,vector):

	diff = True
	for i in mat1:
		if list(i) == list(vector):
			diff = False
	return diff

def FracToCart(a,b,c,alpha,beta,gamma):
	# Calculate transformation matrix to transform fractional coordinates
	# to cartesian

	# Degrees to radians

	alpha = alpha*m.pi/180
	beta = beta*m.pi/180
	gamma = gamma*m.pi/180

	cosa = m.cos(alpha)
	cosb = m.cos(beta)
	cosg = m.cos(gamma)

	cos2a = cosa**2
	cos2b = cosb**2
	cos2g = cosg**2

	sing = m.sin(gamma)

	vol = m.sqrt(1-cos2a-cos2b-cos2g + (2*cosa*cosb*cosg))

	mat = np.ones((3,3))

	mat[0,0] = a
	mat[0,1] = b*cosg
	mat[0,2] = c*cosb

	mat[1,0] = 0.0
	mat[1,1] = b*sing
	mat[1,2] = c*(cosa - (cosb*cosg))/sing

	mat[2,0] = 0.0
	mat[2,1] = 0.0
	mat[2,2] = c*vol/sing

	mat = CleanMatElements(mat)

	# Transpose transformation matrix, so the vectors are aligned in
	# columns, not rows, i.e
	# mat = [a1,a2,a3
	#        b1,b2,b3
	#        c1,c2,c3]
	mat = mat.T
	
	return mat

def CleanMatElements(mat):

	# Clean numeric mess in numpy array, by replacing really small 
	# elements with 0
	# For instance:
	# input
	# [  0.00000000e+00   5.43053000e+00   3.32524059e-16] 
	# output
	# [  0.0   5.43053   0.0]

	small =  abs(mat) < 1e-12
	mat[small] = 0.0

	return mat

def ReadSym(string):
	# Reads the _symmetry_equiv_pos_as_xyz operation in CIF 
	# Assumes that general format is:
	# 'x+y+z+int1/int2'
	# Example: for 'x-y+1/2' it can be written in form of lists:
	# operation=[1,1,0,0.5]  x=1, y=1, z=0, frac=0.5
	# signs=[1,-1,1,1]   +x -y +z +fraction

	operation = [0,0,0,0]
	# operation[0]  -  1 if there is x symmetry operation, 0 otherwise
	# operation[1]  -  1 if there is y symmetry operation, 0 otherwise
	# operation[2]  -  1 if there is z symmetry operation, 0 otherwise
	# operation[3]  -  1 if there is fraction, 0 otherwise
	
	signs = [1,1,1,1] 
	# signs[0] - 'x' sign, -1 means -x
	# signs[1] - 'y' sign, -1 means -y
	# signs[2] - 'z' sign, -1 means -z
	# signs[3] - fraction sign, -1 means -fraction
	
	# Find x,y,z 
	idx = string.find("x")
	if idx > -1: 
		operation[0] = 1
		# Find sign before x
		if idx != 0: #check in the case fraction is first element
			if string[idx-1] == '-':
				signs[0] = -1

	idx = string.find("y")
	if idx > -1: 
		operation[1] = 1
		# Find sign before x
		if idx != 0: #check in the case fraction is first element
			if string[idx-1] == '-':
				signs[1] = -1

	idx = string.find("z")
	if idx > -1: 
		operation[2] = 1
		# Find sign before x
		if idx != 0: #check in the case fraction is first element
			if string[idx-1] == '-':
				signs[2] = -1


	# Find fraction
	idx = string.find("/")
	if idx > -1:
		f1idx = idx-1
		f2idx = idx+1
		frac = float(string[f1idx])/float(string[f2idx])
		operation[3] = frac
		
		# Find sign before fraction
		if f1idx != 0: #check in the case fraction is first element
			if string[f1idx-1] == '-':
				signs[3] = -1

	return signs,operation

def ReadSymFrac(string):
	# This is specyfic routine for CIF file, to transform symmetry 
	# equivalent position to actual true numbers.
	# Works for the cases when the symmetry operation is of the forms:
	# "x", "-x", "x+1/2", "-z+3/4"...

	# There can be only two operation op1*x+op2*ratio, where 
	# op1 and op2 = "+" or "-". Define op list as two "+" operations
	# 
	# Input:
	# - string - symmetry operation as from CIF file, eg. "x+1/2"
	# Output:
	# - coordid - number for coordinate in symmetry operation: 
	#             'x'=0, 'y'=1, 'z'=2
	# - op - signs in symmetry operation, e.g:
	#        "x+1/2" : op=[1,1]
	#        "-x+1/2": op=[-1,1]
	# 	 "-x-1/2": op=[-1,-1] 
	#        "-x"    : op=[-1,1] ....
	# - digits - the digits in the fraction, e.g.:
	#        "1/4" : digits=[1,4]
	#
	# Author: Jakub Kaminski, UCLA, 04/2013


	op = [1,1]

	# There can be two digits defining ratio, e.g. 1/4. 
	# Define digits list with two "0,1" for the start

	digits = [0,1]

	# Find what is the coordiante
	coorindex = string.find("x")
	coordid = 0
	if coorindex == -1:
		coorindex = string.find("y")
		coordid = 1
		if coorindex == -1:
			coorindex = string.find("z")
			coordid = 2
	
	# Check if the coordinate is not negative 
	# and mark it by -1 in op variable
	if coorindex > 0: # just in case to be sure we are not accessing beyond
		          # string length
		optype = string[coorindex-1]
		if optype == "-":
			op[0] = -1
	
	# Check if we add or substract ratio.
	# If we add, op[1]=1, if we substract op[1]=-1
	if coorindex <len(string): # just in case to be sure we are not accesing
		                   # beyond string length 
		optype = string[coorindex+1]
		if optype == "-":
			op[1] = -1
	

	# Now read the ratio
	# Find digit by digit and save them to digit index
	digitindex = 0
	for i in string:
		if i.isdigit():
			digits[digitindex] = float(i)
			digitindex = digitindex + 1
	
	return coordid, op, digits


def CheckParentheses(input):
	# In some CIF file the cell dimensions and anglesare given with 
	# standard deviation in final digits. It is usally the last 
	# number given in parentheses. 

	# CheckParentheses checks is string read from CIF has paretheses and 
	# disregards them.
	#
	# Variables:
	# input is as string
	# Returns float
	#
	# Author: Jakub Kaminski, UCLA, 04/2013

	i = input.find("(")

	if i> -1: # found
		input = input[:i]

	return float(input)

def extractSpaceGroup(symbol):
	# Extract spacegroup symbol from CIF file line
	symbol = symbol.strip("_symmetry_space_group_name_H-M")
	symbol = symbol.strip() # strip any leading spaces
	symbol = symbol.strip("'") # strip quote signs
	symbol = symbol.strip('"') # strip double quotes (just in case)

        # In case symbol contains "_", i.e "I4_1/amd"
        if symbol.find("_") > 0:
                tmp = ""
                for i in symbol:
                        if i != "_":
                                tmp += i
                symbol = tmp

	# In the case the symbol is given with white space, i.e. "F M -3 M"
	symbol = symbol.split()
	symbolOut = ""
	for i in symbol:
		symbolOut += i
	symbolOut = symbolOut.capitalize()

	# In ISCD "S" or ":S" extension on the end of the symbol is equivalent
	# to :1 notation for origin choice for space groups with two origins,
	# .i.e "F d -3 m S"
	# Our notoations follows :1 :2 extension for different origins and 
	# :R :H for different axes, so convert everything to this notation

	if symbolOut[-1] == "s": 
		symbolOut = symbolOut[:-1] +":1"
	
	if symbolOut[-2:] == ":s": 
		symbolOut = symbolOut[:-2] +":1"
	
	# Intel is using follwing convetion in their MadeA software to
	# label group with different origins
	# FD-3MO1 - origin 1; FD-3MO2 - origin 2

	if symbolOut[-2:] == "o1":
		symbolOut = symbolOut[:-2] +":1"

	if symbolOut[-2:] == "o2":
		symbolOut = symbolOut[:-2] +":2"

	# TODO: What are other convetions?

	return symbolOut

def unique2(a,labels):
	# Remove duplicate elements from array keeping it order
	# Not the fastest but will do for now
	# 3.5s to sort 13122 elements array

	dups = set()

	newlist = []
	newlabels = []
	ii = 0

	for elem in a:
		if str(elem) not in dups:
			newlist.append(elem)
			newlabels.append(labels[ii])
			dups.add(str(elem))
		ii += 1

	return np.array(newlist), newlabels

def fneigh(point,setp):

        # Find of nearest neighbours of the "point" in the "setp" of points
        # 
        # Return array of the shape, where
        # - 1st column is the distance of "point" to point p in "setp"
        # - 2nd column is the index of point p in "setp"

        neigh = np.zeros((len(setp),2))
        idx = 0
        for p in setp:
                d = np.linalg.norm(point-p)
                neigh[idx,0] = d
                neigh[idx,1] = idx
                idx += 1

        # Sort the distances according to the distance
        neigh = neigh[np.argsort(neigh[:,0])]

        neigh = CleanMatElements(neigh)

        return neigh

def calcRuntime(tStart,tStop):
	# Converts time between tStop and tStart (given in seconds)
	# to minues:seconds
	# 
	# Outputs: list where first element is number of minutes,
	#          second element number of seconds

	delta = tStop - tStart

	minutes = int(m.floor(delta/60))
	seconds = delta%60
	
	return [minutes,seconds]

class Surface:
	def __init__(self,transM,positions,atoms,atomTyp,Midx):
		self.A = transM[0] # Vector A of the unit cell
		self.B = transM[1] # Vector B of the unit cell
		self.C = transM[2] # Vector C of the unit cell
		self.positions = positions # coordinates of material
		self.atoms = atoms # list with atom labels
		self.Midx = Midx # Miller indices of the plane
		#Define plane vectors and normal
		self.u = np.array((0,0,0))
		self.v = np.array((0,0,0))
		self.n = np.array((0,0,0))
		self.origin = np.array((0,0,0)) 
		self.a = np.array((0,0,0)) # primitive vector
		self.b = np.array((0,0,0)) # primitive vector

		self.unitcell = self.positions.copy() # save unit cell 
		self.unitcellatoms = atoms # save unit cell atom labels
		self.atomTyp = atomTyp
		self.primarea = 0 # area of the primitive cell
		self.phiab = 0 # angle between primtive vectors
		self.norma = 0 # norm of primitive vector a
		self.normb = 0 # norm of primitive vector b

		self.planepos = np.array([[]])#empty array for surface positions
		self.planeatms = [] # empty list for plane atom labels

		self.planeposblk = np.array([[]])# empty array for surface 
			         		 # positions and atoms below
		self.planeatmsblk = [] # empty list for plane atom labels
		                       # and atome below
		self.positionsSmall = np.array([[]]) # center of bulk
		self.positionsSmall = [] # center of bulk atoms
		self.avecs = np.array((0,0,0)) # anticpatory vectors
		self.nneigh = 0 # number of nearst neigbours
		self.exists = False # flag to check if substrate exists

	def __removeperiodic3D(self,pos,vec1,vec2,vec3):

		# Find the atoms in the superlattice that are
		# translations of the other atoms 

		# Rude and not effecient solution 
		# Try more through array interception

		r = range(2)
		uniqlist=[]
		poscheck = []
		for x in r:
			for y in r:
				for z in r:
					if x != 0 or y!= 0 or z!=0:
						if poscheck == []:
							poscheck = pos+x*vec1+y*vec2+z*vec3
						else:
							tmp = pos + x*vec1 + y*vec2+z*vec3
							poscheck = np.vstack\
							          ([poscheck,tmp])

		# Find indicies of unique elements in pos array that 
		# are not there is poscheck array
		ii = 0
		for i in pos:
			uq = True
			for j in poscheck:
				# Smaller accuracy required for z-axis
				# so round it to 5
				# Gives problems otherwise
				if round(i[0],8) == round(j[0],8) and \
				   round(i[1],8) == round(j[1],8) and \
				   round(i[2],8) == round(j[2],8):
					   uq = False
			if uq:
				if ii not in uniqlist:
					uniqlist.append(ii)
			ii += 1
		return uniqlist


	def bulkNEW(self,ncells):
		# 
		# OLD ROUTINE TO CONSTRUCT BULK MATERIAL USING 
		# __UNIQUE ROUTINE TO FIND DUPLICATE ELEMENTS.
		# THE PROBLEM WAS THAT IT DID NOT KEEP THE ORDER 
		# OF THE POSITION ARRAY
		# Consruct bulk surface from "ncells x unit cell"
		#
		# The results is saved in the positions variable

		r = range(-ncells,ncells+1)
		
		# Initial cooridnates and labels 
		posout = self.unitcell.copy()

		# Make unit cell periodic
		perIDX = self.__removeperiodic3D(self.unitcell,self.A,self.B,self.C)
		unitcellPer = self.unitcell.copy()[perIDX]
		unitcellatomsPer = self.unitcellatoms.copy()[perIDX]

                nElems = len(r)
                nAtoms = len(unitcellPer)
                nElems = nElems * nElems * nElems * nAtoms
                posout = np.zeros((nElems,3))
                newlabels = np.zeros((nElems))
                posout[0:nAtoms] = unitcellPer
                newlabels[0:nAtoms] = unitcellatomsPer

                if ncells> 2:
                        middle = len(r)/2
                        closeR = r[middle-2:middle+3]
                else:
                        closeR = r
                # Generate the middle of the bulk first. This is usefull in the later parts of the code to 
                # limit the size of the bulk to small number of atoms, for instance to look for nearest 
                # neighbors 
                i = nAtoms
                for x in closeR:
                        for y in closeR:
                                for z in closeR:
                                        if x != 0 or y!= 0 or z != 0:
                                                posout[i:nAtoms+i] = unitcellPer \
                                                           + x*self.A \
                                                           + y*self.B \
                                                           + z*self.C
                                                newlabels[i:nAtoms+i] = unitcellatomsPer
                                                i += nAtoms

                idxMiddle = i
                # Continue with the rest of bulk
                if ncells >2:
                        for x in r:
                                for y in r:
                                        for z in r:
                                                if x not in closeR or y not in closeR or z not in closeR:
                                                        posout[i:nAtoms+i] = unitcellPer \
                                                                   + x*self.A \
                                                                   + y*self.B \
                                                                   + z*self.C
                                                        newlabels[i:nAtoms+i] = unitcellatomsPer
                                                        i += nAtoms


		#ii = 0
                #for i in posout:
                #       #print newlabels[ii],i[0],i[1],i[2]
                #       print self.atomTyp[newlabels[ii]],i[0],i[1],i[2]
                #       ii+=1
		self.positions = posout
		self.atoms = newlabels
                self.positionsSmall = posout[0:idxMiddle]
                self.atomsSmall = newlabels[0:idxMiddle]

	def bulk(self,ncells):
		# 
		# OLD ROUTINE TO CONSTRUCT BULK MATERIAL USING 
		# __UNIQUE ROUTINE TO FIND DUPLICATE ELEMENTS.
		# THE PROBLEM WAS THAT IT DID NOT KEEP THE ORDER 
		# OF THE POSITION ARRAY
		# Consruct bulk surface from "ncells x unit cell"
		#
		# The results is saved in the positions variable

		r = range(-ncells,ncells+1)
		
		# Initial cooridnates and labels 
		posout = self.unitcell.copy()
		newlabels = self.__addatoms([])

		for x in r:
			for y in r:
				for z in r:
					if x != 0 or y!= 0 or z != 0:
						newpos = self.unitcell \
						           + x*self.A \
							   + y*self.B \
						           + z*self.C
						posout = \
						      np.vstack([posout,newpos])
						newlabels = \
						      self.__addatoms(newlabels)


		# Remove all duplicates

#		self.positions, self.atoms = self.__unique(posout,newlabels)
		self.positions,self.atoms = unique2(posout,newlabels)

	def bulkNAIVE(self,ncells):
		#
		# Consruct bulk surface from "ncells x unit cell"
		#
		# The results is saved in the positions variable

		r = range(-ncells,ncells+1)
		
		# Initial cooridnates and labels 
		posout = self.unitcell.copy()
		newlabels = self.__addatoms([])


		for x in r:
			for y in r:
				for z in r:
					if x != 0 or y!= 0 or z != 0:
						newpos = self.unitcell \
						           + x*self.A \
							   + y*self.B \
						           + z*self.C
						posout,newlabels = \
							self.__uniqueNAIVE\
							 (newpos,posout,\
							 newlabels)

		self.positions = posout
		self.atoms = newlabels
		
	def __addatoms(self, newlabels):
		# Add unit cell atoms labels to atom label list
		for i in range(len(self.unitcellatoms)):
			newlabels.append(self.unitcellatoms[i])
		return newlabels

	def __addatomsNEW(self, nl,base):
		# Add unit cell atoms labels to atom label list
		for i in range(len(base)):
			nl.append(base[i])
		return nl

	def __uniqueNAIVE(self,newmat,oldmat,labels):

		# VERY NAIVE WAY OF FINDING DUPLICATE ROWS IN ARRAY
		# MEMORY AND CPU TIME INEFFECIENT
		# HAS THE ADVANTAGE OF KEEPING THE ORDER
		# WILL WORK FOR NOW, BUT FIND BETTER WAY ASAP
		# 42s to sort 13122 elements array

		out = oldmat.copy()
		#duplicate = False
		ii = 0
		count= 0
		for elem1 in newmat:
			duplicate = False
			for elem2 in oldmat:
				if (elem1[0] == elem2[0]) and \
				   (elem1[1] == elem2[1]) and \
				   (elem1[2] == elem2[2]):
					duplicate = True
					break # break the for loop

			if not duplicate:
				count += 1
				out = np.vstack([out,elem1])
				labels.append(self.unitcellatoms[ii])
			ii += 1
			

		return out,labels

	def __unique(self,a,labels):

		# Check the numpy array with coordiantes for duplicate entries
		# Remove duplicates, remove lables correspoding to duplicated 
		# coordinates. 
		# The returned coordiantes are sorted in different order than 
		# orginals
		#
		# Subroutine take from:
		# http://stackoverflow.com/questions/8560440/removing-duplicate-columns-and-rows-from-a-numpy-2d-array

		# 0.5s to sort 13122 elements array

		order = np.lexsort(a.T)
		a = a[order]
		diff = np.diff(a, axis=0)
		ui = np.ones(len(a), 'bool')
		ui[1:] = (diff != 0).any(axis=1) 

		ii = 0
		newlabels = []
		for idx in order:
			if ui[ii] != False:
				newlabels.append(labels[idx])
			ii += 1

		return a[ui], newlabels

	def construct(self):
		h = self.Midx[0]
		k = self.Midx[1]
		l = self.Midx[2]

		if h != 0 and k == 0 and l == 0:
			# (100) surface
			self.u = self.B.copy()
			self.v = self.C.copy()

			# We take abs(h) as (-100) will be same as (100)
			self.origin = (1.0/abs(h)) * self.A

		elif h == 0 and k != 0 and l == 0:
			# (010) surfae
			self.u = self.A.copy()
			self.v = self.C.copy()
			
			# We take abs(k) as (0-10) will be same as (010)
			self.origin = (1.0/abs(k)) * self.B

		elif h == 0 and k == 0 and l != 0:
			# (001) surface
			self.u = self.A.copy()
			self.v = self.B.copy()

			# We take abs(l) as (00-1) will be same as (001)
			self.origin = (1.0/abs(l)) * self.C

		elif h != 0 and k != 0 and l == 0:
			# (hk0) surface

			# By symmetry (-h,-k,0) is equvalent to (h,k,0)
			if h < 0 and k < 0:
				h *= -1
				k *= -1
			
			self.u = (1./h)*self.A - (1./k)*self.B
			self.v = self.C.copy()

			if h > 0 and k > 0: #110
#				self.origin = (1./k)*self.B
				self.origin = (1./k)*self.B

			if h < 0 and k > 0: #(-110)
				self.origin = (1./k)*self.B + self.A

			if h > 0 and k <-1: #(1-10)
				self.origin = (1./abs(k))*self.B


		elif h != 0 and k == 0 and l != 0:
			# (h0l) surface

			# By symmetry (-h,-k,0) is equvalent to (h,k,0)
			if h < 0 and l < 0:
				h *= -1
				l *= -1

			self.u = (1./h)*self.A - (1./l)*self.C
			self.v = self.B.copy()

			if h > 0 and l > 0: #(101)
				self.origin = (1./l)*self.C

			if h < 0 and l > 0: #(-101)
				self.origin = (1./l)*self.C + self.A

			if  h > 0 and l < -1:
				self.origin = (1./abs(l))*self.C


		elif h == 0 and k != 0 and l!= 0:
			# (0kl) surface

			# By symmetry (-h,-k,0) is equvalent to (h,k,0)
			if k < 0 and l < 0:
				k *= -1
				l *= -1

			self.u = (1./k)*self.B - (1./l)*self.C
			self.v = self.A.copy()

			if k > 0 and l > 0: #(011)
				self.origin = (1./l)*self.C

			if k < 0 and l > 0: #(0-11)
				self.origin = (1./abs(l))*self.C + self.B

			if k > 0 and l < -1: #(01-1)
				self.origin = (1./abs(l))*self.C


		elif h != 0 and k != 0 and l != 0:
			# (hkl) surface

#			self.u = (1.0/gcd_hk) * (k*self.A - h*self.B)
#			self.v = (1.0/gcd_hl) * (l*self.A - h*self.C)

			# The equivalent planes
			# (-1-1-1)=(111)
			# (1-1-1)=(-111)
			# (-1-11)=(11-1)
			# (-11-1)=(1-11)
			# Find all equvalent planes
			counter = 0
			for idx in self.Midx:
				if idx < 0:
					counter += 1

			if counter >= 2:
				h *= -1
				k *= -1
				l *= -1

			gcd_hk = fractions.gcd(abs(h),abs(k))
			gcd_hl = fractions.gcd(abs(h),abs(l))

#			self.u = (1.0/gcd_hk) * (h*self.A - k*self.B)
#			self.v = (1.0/gcd_hl) * (h*self.A - l*self.C)


			self.u =  (1./h)*self.A - (1./k)*self.B
			self.v =  (1./l)*self.C - (1./k)*self.B

			if h > 0 and k > 0 and l > 0:
				self.origin = (1./k)*self.B

			if h < 0 and k > 0 and l > 0:
				self.origin = self.A + (1./k)*self.B

			if h > 0 and k < -1 and l > 0:
				self.origin = (1./abs(k))*self.B

			if h > 0 and k > 0 and l < 0:
				self.origin = self.C + (1./k)*self.B


	def plane(self):
		
		# Cut the plane from the postitions
		
		# Define normal to the plane
		self.n = normal(self.u,self.v)

		# Shift the orgin of the unit cell so it 
		# corresponds to the plane vectors
		self.positions = self.positions - self.origin

		# Rotate the surface to xy plane
		# Algorithm:
		# 1) Find angle of normal of the plane to the z axis
		# 2) Define vecor of rotation as cross product between
		#    normal and z axis
		# 3) Define rotation matrix for this angle and vector
		# 4) Rotate the plane

		# Find the angle
		normN = np.linalg.norm(self.n) # norm of plane normal

		z = np.array((0,0,1))
		normZ = 1 # norm of Z
		
		cosphi = np.dot(self.n,z)/(normN*normZ)
		phi = m.acos(cosphi)

		phi = phi * 180/m.pi
		
		# Create diagonal rotation matrix 
		R = np.diag((1.0,1.0,1.0))

		if phi != 0: # angle 0 gives problems with NaN in R matrix,
			     # but we dont need it anyways in this case
			phi = 180 - phi # Rotation is counterclockwise, 
			                #so define angle as 180 - phi
			phi = phi * m.pi/180

			cosphi = m.cos(phi)
			sinphi = m.sin(phi)

			rv = np.cross(self.n, z) # rotation vector
			normrv = np.linalg.norm(rv)
			rv = rv/normrv
			
			# now define rotation matrix
#			R = np.zeros((3,3))
			R[0,0] = cosphi + (rv[0]**2)*(1-cosphi)
			R[0,1] = rv[0]*rv[1]*(1-cosphi) - rv[2]*sinphi
			R[0,2] = rv[0]*rv[2]*(1-cosphi) + rv[1]*sinphi
	
			R[1,0] = rv[1]*rv[0]*(1-cosphi) + rv[2]*sinphi
			R[1,1] = cosphi + (rv[1]**2)*(1-cosphi)
			R[1,2] = rv[1]*rv[2]*(1-cosphi)-rv[0]*sinphi
	
			R[2,0] = rv[2]*rv[0]*(1-cosphi) - rv[1]*sinphi
			R[2,1] = rv[2]*rv[1]*(1-cosphi) + rv[0]*sinphi
			R[2,2] = cosphi + (rv[2]**2)*(1-cosphi)

		# Rotate normal
		self.n = np.dot(self.n,R)
		self.n = CleanMatElements(self.n)

		# Rotate plane vectors
		self.u = np.dot(self.u,R)
		self.v = np.dot(self.v,R)
		self.u = CleanMatElements(self.u)
		self.v = CleanMatElements(self.v)
		# Rotate positions of atoms
		posrot = np.dot(self.positions,R)
		posrot = CleanMatElements(posrot)

		# list of indexes of atom belonging to the plane
		idxlist = []
		# list of indexes of atom belonging to the plane and below
		idxlistblk = []

		# Point is on the plane if its vector is orthogonal 
		# to the normal of the plane.

		# TODO: the implementation via looking dot product 
		# with normal is from 
		# old version of the code. Now its redundant, and can 
		# be simplified by specyfing that atoms in plane are the
		# ones with index z = 0 , above the plane are the ones
		# with z > 0 and below with z < 0 

		idxlistBOOL = posrot[:,2] == 0
		idxlist = np.where(idxlistBOOL == 1)[0]
		idxlistblkBOOL = posrot[:,2] < 0
		idxlistblk = np.where(idxlistblkBOOL == 1)[0]
		
		if idxlist != []: self.exists = True

		if self.exists:

			# Create anticipatory vectors
			# Find all atoms lying above the plane
			idxabove = posrot[:,2]>0
			posabove = posrot[idxabove]

			# Labels of atoms on surface only
			self.planeatms = self.atoms[idxlist]
	
			uqatoms,uqidx,uqidxbulk = self.__finduqatoms(self.planeatms,\
					                   posrot,idxlist)
				
			# Find anticipatory vectors for each unique atom on the plane
			avecslist = []
			for idx in uqidxbulk:
				avecsitem = self.__anticipatoryvecs\
   					     (posrot[idx],posabove)
				avecslist.append(avecsitem)

			self.avecs = avecslist[0]
			# Create a dictionary avecall to hold all anticipatory
			# vectors assiciated to given atom type. This will be 
			# usefull for scroing function
			self.avecsall = {}
			for i in range(len(uqatoms)):
				lab = uqatoms[i]
				ii = avecslist[i]
				self.avecsall[lab]=ii
	
#			# Find nearest unique atoms in the whole bulk 
#			# This is needed to find nearest neigbours of each atom type
#			uqatoms,uqidx,uqidxbulk = self.__finduqatoms(self.atoms,\
#					                   posrot,\
#							   idxlist=range(len(posrot)))
                        uqatoms,uqidx,uqidxbulk = self.__finduqatoms(self.atomsSmall,\
                                                           self.positionsSmall,\
                                                           idxlist=range(len(self.positionsSmall)))

                        neighlist = []
                        for i in uqidxbulk:
                                # Construct array with poistion without atom i
                                postmpa = self.positionsSmall[:i]
                                postmpb = self.positionsSmall[i+1:]
                                postmp  = np.concatenate((postmpa,postmpb))
                                nneighat = self.__anticipatoryvecs(self.positionsSmall[i],\
                                                    postmp,neighonly=True)
                                neighlist.append(nneighat)

			# Create dictionary that as a key will have label 
			# of unique atom and as value, number of its nearest 
			# neighbors
			self.uqneigh = {}
			for i in range(len(uqatoms)):
				lab = uqatoms[i]
				ii = neighlist[i]
				self.uqneigh[lab]=ii
			self.nneigh = neighlist[0]
			
			# Revert to orginal origin of atom coordinates
			self.positions = self.positions + self.origin
	
			# Create surface coordinates
			self.planepos = posrot[idxlist]
			
			# Create surace coordiantes including atoms below it
			# Make sure that surface atoms are first in the array
			#idxlistblk = idxlist + idxlistblk
			idxlistblk = np.concatenate((idxlist,idxlistblk))
			self.planeposblk = posrot[idxlistblk]

			# Labels of atoms on surface and below
			self.planeatmsblk = self.atoms[idxlistblk]
	
			# Translate plane so coordinate z=0	
#			if abs(self.planepos[:,2]).all > 0:
#				self.planepos[:,2] -= self.planepos[:,2]

			# end if self.exists

	def __finduqatoms(self,labels,pos,idxlist):
		# Routine to find unique types of atoms from set of
		# atom lables
		# Return atom labels, and index of the atom on the plane
		# and index of atoms in the bulk structure

		uql = [] # unique labels
		uqi = [] # index of representative atom on surface
		uqibulk = [] # index of representative atom in bulk structure

		postmp = pos[idxlist]	
		for i in range(len(postmp)):
			atom = labels[i]
			if atom not in uql:
				uql.append(atom)
				uqi.append(i)

		for i in uqi:
			uqibulk.append(idxlist[i])

		return uql,uqi,uqibulk

	def __anticipatoryvecs(self,atom,bulk,neighonly=False):
		# Routine to find anticipatory vectors for atom
		# by looking through it nearest neighbours 

		# Find nearest neighbours of the atom
		# Construct array with the distance and indices of atoms in 
		# bulk
		# Introcude threshold to find number neartest neighbours
		# There are cases where atoms has N nearest neighbours, but
		# but they have only slighlty different distances, for instance
		# in the case of a-quarts-SiO2, Si has four NN, two has distance
		# 1.62A and two 1.63A. Use threshold to find them
		Nthresh = 0.1

		neigh = np.zeros((len(bulk),2))
		idx = 0
		for atm in bulk:
			d = np.linalg.norm(atom-atm)
			neigh[idx,0]=d
			neigh[idx,1]=idx
			idx += 1

		# Sort the distances according to the distance
		neigh = neigh[np.argsort(neigh[:,0])]

		# To avoid small differences in floating point comparisons
		# substract first elements from all others and round the 
		# result to 12 decimal place
		shift = neigh[0,0].copy()
#		neigh[:,0] -= neigh[0,0] 
		neigh[:,0] -= shift

		neigh = CleanMatElements(neigh)
		
		# Select only atoms with the shortest distance
		#idx = neigh[:,0] == neigh[0,0]
		#idx = neigh[:,0] == 0.0
		idx = neigh[:,0] <= Nthresh
		neigh = neigh[idx]

		if neighonly:
			return len(neigh)

		# Find the coordinates of nearest neighbours in bulk
		idxs = neigh[:,1].astype(int)
		avecs = bulk[idxs]
		
		# Shift avecs to origin of coordinate system
		avecs = avecs - atom
		avecs = CleanMatElements(avecs)

		# If there is more than one anitcipatory vector,
		# see if they are quivalent. If yes, remove the duplicates
		# The vectors are  assumed to be equivalent if the angle
		# between vectors and z-axis is the same
		
#		if len(avecs) > 1:
#			z = np.array((0,0,1))
#			nz = 1
#			ang = []
#			idx = []
#			ii = 0
#			for v in avecs:
#				nv = np.linalg.norm(v)
#				cosphi = np.dot(v,z)/(nv*nz)
#				if cosphi not in ang:
#					ang.append(cosphi)
#					idx.append(ii)
#				ii += 1
#			avecs = avecs[idx]

		return avecs

	def initpvec(self):
		# Initialize primitive vecotrs
		# Algorithm:
		# Take first atom in the surface, calculate its nearest 
		# neighbours, and define primitive vectors as the 
		# vectors pointing to the nearest two atoms. 

		# Shift coordinates so that 1st atom is in origin
		orgsave = self.planepos[0]
		self.planepos -= self.planepos[0]

		# Find all the nearest neighbours from atom 0
		distmat = np.sum(np.abs(self.planepos)**2,axis=1)**(1./2)
		# Sort distmat. 
		idx = distmat.argsort()

		# Find two non-linear vectors. They will be initial 
		# primitive vectors
		# It is most likely to be fist two smallest vectors,
		# but check for linearity just in case

		linear = True
		i = 1
		while i < len(distmat):
			j = i + 1
			while j < len(distmat):
				self.a = self.planepos[idx[i]]
				self.b = self.planepos[idx[j]]

				# In the case of crystals containig different
				# atoms, the primitive vecors need to be defined
				# between atoms of the same type.
				if self.planeatms[idx[i]] == self.planeatms[0]\
				and self.planeatms[idx[j]] == self.planeatms[0]:

					check = np.cross(self.a,self.b)
					# Clean numeric noise
					check = CleanMatElements(check)
					#TMP
					self.norma=np.linalg.norm(self.a)
					self.normb=np.linalg.norm(self.b)
					print "NORM PRIM ",self.norma,self.normb
					cosphi = np.dot(self.a,self.b)/(self.norma*self.normb)
					if np.linalg.norm(check) != 0:
						linear = False

				if not linear: break
				j += 1
			i += 1
			if not linear: break

                if linear:
                        print "Only linear primitive vectors found"
                        print "Something is wrong. Exiting."
                        exit()

		# Reduce a and b
		self.a, self.b = reduction(self.a, self.b)

#		print "REDUCED VECTORS"
#		print "A", self.a, np.linalg.norm(self.a)
#		print "B", self.b, np.linalg.norm(self.b)
		
		# Shift coordinates back to orginals positions 
		self.planepos += orgsave

	def initpvecNEW(self):
		# Initialize primitive vecotrs
		# Algorithm:
		# Take first atom in the surface, calculate its nearest 
		# neighbours, and define primitive vectors as the 
		# vectors pointing to the nearest two atoms. 

		# Algortihm used:
		# Check two conditions:
		# 1) If lattice constant for given plane can be expresed 
		#    in term of intiger multiplications of primitive vector:
		#        n*|a| = |u| , n=1,2,3,....
		# 2) If length of the scalar projection of the primitive 
		#    vector on the lattice vector is equal to 0.5*(|u|**2)
		#    This can be shown from the properties of dot product:
		#   dot(a,u) == 0.5*(|u|**2) when a_u = |a|*cos(phi) == 0.5*|u|
		#
		#   If any of those conditions is met, the pair of vectors 
		#   are the primitive vectors for this lattice

		# Primitive vectors not found yet
		self.exists = False
		
		# Shift coordinates so that 1st atom is in origin
		orgsave = self.planepos[0].copy()
		self.planepos -= orgsave

		# Find all the nearest neighbours from atom 0
		distmat = np.sum(np.abs(self.planepos)**2,axis=1)**(1./2)
		# Sort distmat. 
		idx = distmat.argsort()

		# Find two non-linear vectors. They will be initial 
		# primitive vectors
		# It is most likely to be fist two smallest vectors,
		# but check for linearity just in case

		nu = np.linalg.norm(self.u)
		nv = np.linalg.norm(self.v)

		primitive = False

		# Search for 1st primitive vector
		i = 1
		while i < len(distmat):
			# In the case of crystals containig different
			# atoms, the primitive vecors need to be defined
			# between atoms of the same type.
			if self.planeatms[idx[i]] == self.planeatms[0]:
				self.a = self.planepos[idx[i]].copy()
	
				na = np.linalg.norm(self.a)
	
				# 1st condtition
				if nu >= na:
					ma = nu%na
				else:
					ma = na%nu
	
#				print nu,na,ma
	
				if round(ma,6) == 0.0: primitive = True
	
				# 2nd condition
				dau = np.dot(self.a,self.u)
	
				if round(abs(dau),5) == round((nu**2)/2,5): primitive = True

				if primitive: break
			i += 1

		# Search for 2nd primitive vector
		linear = True
		primitive = False
		i = 1
		while i < len(distmat):
			# In the case of crystals containig different
			# atoms, the primitive vecors need to be defined
			# between atoms of the same type.
			if self.planeatms[idx[i]] == self.planeatms[0]:
				self.b = self.planepos[idx[i]].copy()
	
				nb = np.linalg.norm(self.b)
	
				# 1st condtition
				if nv >= nb:
					mb = nv%nb
				else:
					mb = nb%nv
	
				if round(mb,6) == 0.0 : primitive = True
	
				# 2nd condition
				dbv = np.dot(self.b,self.v)
				
				if round(abs(dbv),5) == round((nv**2)/2,5): primitive = True
	
				# Check if vector b is not linear with vector a
				if primitive:
					check = np.cross(self.a,self.b)
					# Clean numeric noise
					check = CleanMatElements(check)
					self.norma=np.linalg.norm(self.a)
					self.normb=np.linalg.norm(self.b)
					cosphi = np.dot(self.a,self.b)/\
						       (self.norma*self.normb)
					if np.linalg.norm(check) != 0:
						linear = False
	
				if primitive and (not linear): break

			i += 1
		
		if primitive:
			self.exists = True

			# Reduce a and b
			self.a, self.b = reduction(self.a, self.b)

		# Shift coordinates back to orginals positions 
		self.planepos += orgsave

	def initpvecNEW2(self):
		# Initialize primitive vecotrs
		# Algorithm:
		# Take first atom in the surface, calculate its nearest 
		# neighbours, and define primitive vectors as the 
		# vectors pointing to the nearest two atoms. 

		# Algortihm used:
		# Check two conditions:
		# 1) If lattice constant for given plane can be expresed 
		#    in term of intiger multiplications of primitive vector:
		#        n*|a| = |u| , n=1,2,3,....
		# 2) If length of the scalar projection of the primitive 
		#    vector on the lattice vector is equal to 0.5*(|u|**2)
		#    This can be shown from the properties of dot product:
		#   dot(a,u) == 0.5*(|u|**2) when a_u = |a|*cos(phi) == 0.5*|u|
		#
		#   If any of those conditions is met, the pair of vectors 
		#   are the primitive vectors for this lattice

		# Primitive vectors not found yet
		self.exists = False

		# Shift coordinates so that 1st atom is in origin
		orgsave = self.planepos[0].copy()
		self.planepos -= orgsave

		# Find all the nearest neighbours from atom 0
		distmat = np.sum(np.abs(self.planepos)**2,axis=1)**(1./2)
		# Sort distmat. 
		idx = distmat.argsort()

		# Find two non-linear vectors. They will be initial 
		# primitive vectors
		# It is most likely to be fist two smallest vectors,
		# but check for linearity just in case

		nu = np.linalg.norm(self.u)
		nv = np.linalg.norm(self.v)

		primitive = False

		# Search for 1st primitive vector
		i = 1
		while i < len(distmat):
			# In the case of crystals containig different
			# atoms, the primitive vecors need to be defined
			# between atoms of the same type.
			if self.planeatms[idx[i]] == self.planeatms[0]:
				self.a = self.planepos[idx[i]].copy()
	
				na = np.linalg.norm(self.a)
	
				# 1st condtition
				if nu >= na:
					ma = nu%na
				else:
					ma = na%nu
	
				if round(ma,6) == 0.0: primitive = True
	
		#		# 2nd condition
		#		dau = np.dot(self.a,self.u)
		#		print "2nd condition:",round(abs(dau),5),round((nu**2)/2,5)
		#		if round(abs(dau),5) == round((nu**2)/2,5): primitive = True

				# 3rd condition u and a are colinear
				check = np.cross(self.a,self.u)
				# Clean numeric noise
				check = CleanMatElements(check)
				cosphi = np.dot(self.a,self.u)/\
					       (na*nu)
				if np.linalg.norm(check) == 0:
					primitive = True


				if primitive: break
			i += 1
		# Search for 2nd primitive vector
		linear = True
		primitive = False
		i = 1
		while i < len(distmat):
			# In the case of crystals containig different
			# atoms, the primitive vecors need to be defined
			# between atoms of the same type.
			if self.planeatms[idx[i]] == self.planeatms[0]:
				self.b = self.planepos[idx[i]].copy()
	
				nb = np.linalg.norm(self.b)
	
				# 1st condtition
				if nv >= nb:
					mb = nv%nb
				else:
					mb = nb%nv
	
				if round(mb,6) == 0.0 : primitive = True
	
				# 2nd condition
		#		dbv = np.dot(self.b,self.v)
		#		
		#		if round(abs(dbv),5) == round((nv**2)/2,5): primitive = True
		#
				# 3rd condition u and a are colinear
				check = np.cross(self.b,self.v)
				# Clean numeric noise
				check = CleanMatElements(check)
				cosphi = np.dot(self.b,self.v)/\
					       (nb*nv)
				if np.linalg.norm(check) == 0:
					primitive = True

				# Check if vector b is not linear with vector a
				if primitive:
					check = np.cross(self.a,self.b)
					# Clean numeric noise
					check = CleanMatElements(check)
					self.norma=np.linalg.norm(self.a)
					self.normb=np.linalg.norm(self.b)
					cosphi = np.dot(self.a,self.b)/\
						       (self.norma*self.normb)
					if np.linalg.norm(check) != 0:
						linear = False
	
				if primitive and (not linear): break

			i += 1
		
		if primitive:
			self.exists = True

			# Reduce a and b
			self.a, self.b = reduction(self.a, self.b)

		# Shift coordinates back to orginals positions 
		self.planepos += orgsave


	def primitivecell(self):
		# Calculate the norm of primitive vectors a,b , 
		# angle between them and the area of primitive cell

		self.norma=np.linalg.norm(self.a)
		self.normb=np.linalg.norm(self.b)

		cosphi = np.dot(self.a,self.b)/(self.norma*self.normb)
		cosphi = round(cosphi,12)
		phi = m.acos(cosphi)
		self.phiab = phi * 180/m.pi
		sinphi = m.sin(phi)

		self.primarea = self.norma * self.normb * sinphi

def normal(vec1,vec2):
	# Calculate normal to the plane given by vec1 and vec2
		result = np.cross(vec1,vec2)
		return result

def rotmat2d(phi):
	# Calculate clockwise rotation matrix around z-axis 
	# 
	# phi - input angle, in radians

	mat = np.zeros((2,2))

	mat[0,0] =  m.cos(phi)
	mat[0,1] = -1*m.sin(phi)
	mat[1,0] =  m.sin(phi)
	mat[1,1] =  m.cos(phi)

	return mat


def FindAllDivisors(x):
	# Algorithm taken from
	# http://stackoverflow.com/questions/12421969/finding-all-divisors-of-a-number-optimization
	divList = []
	y = 1
	while y <= m.sqrt(x):
		if x % y == 0:
			divList.append(y)
			# If x is a square of y, don't add it to div list to
			# avoid doubling the same divisor
			if (int(x/y)) != y:
				divList.append(int(x / y))
		y += 1
	divList.sort()
	return divList

def reduction(a,b):
	# Primitive vectors reduction algortihm
	# Refeference:
	# Zur, McGil, J. Appl. Phys. 55, 378 (1984)

	reduced = False

	while not reduced:
		dot = np.dot(a,b)

		if dot < 0:
			b = -1 * b

		norma = np.linalg.norm(a)
		normb = np.linalg.norm(b)

		if norma <= normb:
			ab = a+b
			normab = np.linalg.norm(ab)

			if normb <= normab:
				amb = a - b
				normamb = np.linalg.norm(amb)
				if normb <= normamb:
					reduced = True
				else:
					b = b - a
			else:
				b = b + a
		else: 
			tmp = a
			a = b
			b = tmp
	
	return a,b



def Superlattice(a,b,n):

	# Find all the possible superlattices formed from primitive cell
	# (a,b) multiplied by n
	# a - primitive cell vector a
	# b - primitive cell vector b
	# n - primitive cell multiplicator

	# Based on  Zur, McGil, J. Appl. Phys. 55, 378 (1984) eqs. (2.3)-(2.6)

	# Find the divisors on n

	divlist = FindAllDivisors(n)

	# Construct transformation (2.3) from Zur-McGil
	results = []
	for m in divlist:
		i = n/m # Eq. (2.4)
		for j in range(m):
			tmat = np.zeros((2,2))
			tmat[0,0] = i
			tmat[0,1] = j
			tmat[1,1] = m

			vec = np.array((a[0:2],b[0:2]))

			result = np.dot(tmat,vec)

			# Reduce pair of superlattice vectors
			reda, redb = reduction(result[0], result[1])

			# store the result as vector in 3D (with z = 0)
			tmp1 = np.array((0.0,0.0,0.0))
			tmp2 = np.array((0.0,0.0,0.0))

			tmp1[0:2] = tmp1[0:2] + reda
			tmp2[0:2] = tmp2[0:2] + redb

			pair = np.array((tmp1,tmp2))

			results.append(pair)
	
	return results

class Interface:
	def __init__(self,vecSub,Substrate,nL):

		# Substrate surface and atoms below it
		posSubBlk = Substrate.planeposblk.copy() 
		# Substrate atom labels
		SubAtomsBlk = Substrate.planeatmsblk
		#Number of nearest neigbors of the atom in the bulk of Substrate
		SubNneigh = Substrate.nneigh	

		SubNneighuq = Substrate.uqneigh

		self.vecSub = vecSub.copy() # vectors defining substrate surface

		# Prepare xyz output for calculations
		periodic = True
		nlayers = nL

		self.SSurfPos,self.SSurfAtm,self.vecSubR, =\
		 	 			self.__CreateSurface\
				                (posSubBlk,SubAtomsBlk,\
				                 self.vecSub[0],self.vecSub[1],\
						 nlayers,periodic)

	def __CreateSurface(self,posin,atoms,vec1,vec2,nlayers,\
			    periodic=True):
		# Create surface given reduced superlattice vectors 
		# and plane coordinates coordinates
		# The output surface is rotated so that 
		# vec1 point x direction
		# Return
		# - positions and labels of surface
		# - rotated anticipatory vectors
		# - rotated vectors
		# list of indexes of atom belonging to the plane
		idxlist = []
		planeatms = []

		pos = posin.copy()

		# Limit only to nlayers of atoms
		iidx = abs(pos[:,2]) <= nlayers
		pos = pos[iidx]

		# Generate atom labels for limited positions
		tmplabels = atoms[iidx]
		idx = 0

                # Orient the atoms in such a way that the coorindate origin is 
                # in the middle of the top plane
                # Find positions of only top plane
                idxPlane = abs(pos[:,2]) == 0.0
                posPlane = pos[idxPlane]
                # Find centroid of this plane
                cX = sum(posPlane[:,0])/len(posPlane[:,0])
                cY = sum(posPlane[:,1])/len(posPlane[:,1])
                # The z-variable of centroid is constant
                cZ = posPlane[0,2]
                cXYZ = np.array((cX,cY,cZ))
                nnXYZ = fneigh(cXYZ,posPlane)
                originIdx = int(nnXYZ[0][1])

		# Find the all the atoms that are inside the area 
		# designated by vectors vec1 and vec2
		# We will use barycentric coordinates to do this
		# by defining two triangles with following corners:
		# 1) (0,0,0), vec1, vec2
		# 2) vec1+vec2, vec1, vec2
	
		# Shift the coordiante system so the origin is on first atom
		# We need this for calculations in barycentric coordinates

#		pos -= pos[0] # Does not work for large number of atoms 
		shift = pos[0].copy()
		shift = posPlane[originIdx].copy()
		pos -= shift

		for atom in pos:

			# Lower triangle
			# p1 = [0,0]
			# p2 = vec1
			# p3 = vec2
			p1 = np.array((0.0,0.0))

			alpha1,beta1,gamma1 = self.__barycentric(atom,p1,\
						              vec1,vec2)

			a1 = alpha1 >=0 and alpha1 <= 1
			b1 = beta1  >=0 and beta1  <= 1
			c1 = gamma1 >=0 and gamma1 <= 1

			# Upper triangle
			# p1 = vec1+vec2
			# p2 = vec1
			# p3 = vec2

			p1 = vec1+vec2

			alpha2,beta2,gamma2 = self.__barycentric(atom,p1,\
						              vec1,vec2)


			a2 = alpha2 >= 0 and alpha2 <= 1
			b2 = beta2  >= 0 and beta2  <= 1
			c2 = gamma2 >= 0 and gamma2 <= 1

			if (a1 and b1 and c1) or (a2 and b2 and c2):
				idxlist.append(idx)

			idx += 1

		# Create the surface
		planepos = pos[idxlist]

		# and corresponind atom labels
		planeatms = tmplabels[idxlist]

#		ii=0
#		for i in planepos:
#			print planeatms[ii],i[0],i[1],i[2]
#			ii+=1

		# Rotate the plane in such a way that vec1 is aligned 
		# with x axis
		# Find the angle
		vecx = np.array((1.0,0.0,0.0))
		normvec1 = np.linalg.norm(vec1)
		normx = np.linalg.norm(vecx)
		cosphi = np.dot(vec1,vecx)/(normvec1*normx)
		phi = m.acos(cosphi)
		
		# If vec1 is in 1st and 2nd quater of coordinate system
		# define the angle such that  rotation is clockwise to x axis
		if vec1[1] > 0:
	        	phi =  phi * 180/m.pi
		       	phi = 360 - phi
		        phi = phi * m.pi/180
		
		# Find the rotation matrix
		rotmat = rotmat2d(phi)

		# Rotate plane positions so vector u is aligned with x
		tmp = np.dot(rotmat,planepos[:,:2].T)
		planepos[:,:2] = tmp.T
		
		planepos = CleanMatElements(planepos)

		# Rotate vectors
		vec1r = np.array((0.0,0.0,0.0))
		vec2r = np.array((0.0,0.0,0.0))
		vec1r[:2] = np.dot(rotmat,vec1[:2])
		vec2r[:2] = np.dot(rotmat,vec2[:2])
		vec1r = CleanMatElements(vec1r)
		vec2r = CleanMatElements(vec2r)
		# If vec2 is pointing toward 3th and/or 4th quater of coordinate
		# system, flip the coordinates so they are in 1st and 2st 
		# quater 
		if vec2r[1] <0:
			planepos[:,1] *= -1
			vec2r[1] *= -1

		if periodic: 
		# Find only unique atoms
			uniqueidx = self.__removeperiodic(planepos,vec1r,vec2r)
			planepos = planepos[uniqueidx]
			planeatms = planeatms[uniqueidx]

		# Shift coordinates so they are above 0.0, 0.0, plane
		#minZ = min(planepos[:,2])
		#planepos[:,2] -= minZ
		# The cell is starts and 0.0, 0.0,0.0 and goes below.
		# Flib it so that the bottom surface is at the top, 
		# and other surface stays the same
		planepos[:,2] *= -1

		return planepos, planeatms, [vec1r,vec2r]


	def __barycentric(self,p,p1,p2,p3):

		# Converstion to barycentric coordinates
		# Source:
		# http://en.wikipedia.org/wiki/Barycentric_coordinate_system_%28mathematics%29
		# TODO:
		# BE CAREFUL ABOUT THE ROUNDING!
		alpha = ((p2[1] - p3[1])*(p[0] - p3[0]) +\
			 (p3[0] - p2[0])*(p[1] - p3[1])) /\
			((p2[1] - p3[1])*(p1[0] - p3[0]) + \
			 (p3[0] - p2[0])*(p1[1] - p3[1]))

		beta =  ((p3[1] - p1[1])*(p[0] - p3[0])  +\
			 (p1[0] - p3[0])*(p[1] - p3[1])) /\
			((p2[1] - p3[1])*(p1[0] - p3[0]) + \
			 (p3[0] - p2[0])*(p1[1] - p3[1]))

		# Clean numeric noise - round to 6th place
		alpha = round(alpha,6)
		beta = round(beta,6)

		gamma = 1.0 - alpha - beta
		gamma = round(gamma,6)
		
		return alpha, beta, gamma

	def __removeperiodic(self,pos,vec1,vec2):

		# Find the atoms in the superlattice that are
		# translations of the other atoms 

		# Rude and not effecient solution 
		# Try more through array interception

		r = range(2)
		uniqlist=[]
		poscheck = []
		for x in r:
			for y in r:
				if x != 0 or y!= 0:
					if poscheck == []:
						poscheck = pos+x*vec1+y*vec2
					else:
						tmp = pos + x*vec1 + y*vec2
						poscheck = np.vstack\
						          ([poscheck,tmp])

		# Find indicies of unique elements in pos array that 
		# are not there is poscheck array
		ii = 0
		for i in pos:
			uq = True
			for j in poscheck:
				# Smaller accuracy required for z-axis
				# so round it to 5
				# Gives problems otherwise
				if round(i[0],8) == round(j[0],8) and \
				   round(i[1],8) == round(j[1],8) and \
				   round(i[2],8) == round(j[2],8):
					   uq = False
			if uq:
				if ii not in uniqlist:
					uniqlist.append(ii)
			ii += 1
		return uniqlist

    def checkPolar(self,atom_coor,atomLabels, z_digit=4, dis_tol_rate=0.01, max_compare = float("Inf")):
        # Description of the input variables:
        # - self: not used here, but required by Python
        # - atom_coor: of the atoms, numpy.array(Natoms,3)
        # - atomLabels: numpy.array with integers denoting atom types, 
        #   i.e [0,1,1,0,2,...,natoms]. The order of the atoms is the same 
        #   as in positions array.
        # - atomTypes: dictionary that allows to decode entries in the atomLabels
        #   in terms of real chemical species. 
        #   Example:
        #    atomTypes = {0: 'Ga', 1: 'As', 2: 'H'}
        #    which means that integer 0 corresponds to "Ga", integer 1 to "As" and 2 to "H"
        #   Usage:
        #    find what is the atom type of the 3rd atom in the structure:
        #    atomLabels = [0,1,1,0,2]
        #    atom = atomLabels[2]  # remeberin Python we count from 0, so 3rd atom is 2nd in the structure
        #    type = atomTypes[atom]
        #    In this case atom will be set to "1", and type to "As"
        #
        # - vecX: lattice vector X, numpy.array[x1, x2, x3]
        # - vecY: lattice vector Y, numpy.array[y1, y2, y3]

        #  Start implementation:
        #  Do lots of cool stuff here

        #  Return values are 
        #  polar = {True,False} - is structure polar or not
        #  periodicity = 0  - double number saying what is the periodicity of the strucure in angstrom
        # minDiffRate - double number saying what is the minimal distance rate found during the chekcing procedure (especially usful for understanding the polar structures)
        
        # Setting dummy values now:
        polar = False
        periodicity = float("nan")
        minDiffRate = -float("Inf")

        vecX = self.vecSubR[0,range(2)];
        vecY = self.vecSubR[1,range(2)];

        # first 180 primes
        primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069];

        nAtoms = len(atom_type);
        # represent the atom types by prime numbers
        ele_n = copy.deepcopy(atomLabels); 
        for i in xrange(min(atomLabels), max(atomLabels+1)):
            ele_n = np.where(atomLabels == i, primes[i], ele_n);
        ele_n = outer(atomLabels,atomLabels);
        ele_n[range(nAtoms), range(nAtoms)] = 0;
        
        # build the distance matrix (size nAtoms*nAtoms, entry i,j represent for the distance between i-th atom and j-th atom)
        atom_dist = np.zeros([nAtoms, nAtoms], dtype=float);
        for i in xrange(0, nAtoms):
            for j in xrange(i,nAtoms):
                atom_dist[i,j] = compute_min_dist(atom_coor[i, 0:2] - atom_coor[j, 0:2], vecX, vecY);
                atom_dist[i,j] = atom_dist[i,j]+(atom_coor[i, 2] - atom_coor[j, 2])^2;
                atom_dist[j,i] = atom_dist[i,j]

        # round the z-coordinate
        atom_coor = np.floor(atom_coor[:, 2]*(10^z_digit)) /(10**z_digit);
        # atoms with same cut_z coord are considered on the same "surface"
        # return_inverse=False: the resut is ordered from small to large
        surface_z= np.unique(atom_coor[:, 2], return_inverse=False);  
        surface_n=len(surface_z, 1); # number of surfaces
        init_z = surface_z[-1]; # the first surface to consider
        firstRound = False;
        forDelete = np.zeros([nAtoms, 1], dtype=bool);
        while(surface_n>=2 and not((abs(surface_z[0]-surface_z[-1])<abs(init_z - surface_z[-1])) and not(firstRound))): # the first round non-polar must be found within the upper half of the structure
            if max_compare> (surface_n/2): # take advange of integer division
                max_compare=(surface_n/2);
            polar=False;
            doCompare = True;
    
            if firstRound:
                if abs(init_z -  surface_z[-1]) < perodicity_lower_bound:
                    doCompare = False;

            if doCompare:
                for nSurface in xrange(0,max_compare):
                #each row represents for an atom
                # the indices of atoms on the upper surface & lower surface
                    u_lidx = atom_coor[:, 2]==surface_z[nSurface];
                    u_lidx = u_lidx [~forDelete];
                    l_lidx =atom_coor[:, 2]==surface_z[-nSurface-1];
                    l_lidx = l_lidx[~forDelete];
                
                
                    # if the number of atoms are differnent
                    if sum(u_lidx) != sum(l_lidx):
                        polar=True;
                        break;
                       
                    nAtomPerLayer = sum(u_lidx);
                    data_upper_type= ele_n[u_lidx,:];
                    data_lower_type= ele_n[l_lidx,:];
                    data_upper_dist= atom_dist[u_lidx,:];
                    data_lower_dist= atom_dist[l_lidx,:];
                
                    # for each atom, sort the distance from this atom to the others (small to large)
                    sort1idx_u = np.argsort(data_upper_dist, axis = 1);
                    sort1idx_l = np.argsort(data_lower_dist,axis = 1);
                    for i in xrange(nAtomPerLayer):
                        data_upper_type[i,:] = data_upper_type[i,sort1idx_u[i,:]];
                        data_lower_type[i,:] = data_lower_type[i,sort1idx_l[i,:]];
                        data_upper_dist[i,:] = data_upper_dist[i,sort1idx_u[i,:]];
                        data_lower_dist[i,:] = data_lower_dist[i,sort1idx_l[i,:]];
                    # for each atom, sort the type from this atom to the others (small to large)
                    sort2idx_u = np.argsort(data_upper_type, axis = 1);
                    sort2idx_l = np.argsort(data_lower_type,axis = 1);
                    for i in xrange(nAtomPerLayer):
                        data_upper_type[i,:] = data_upper_type[i,sort2idx_u[i,:]];
                        data_lower_type[i,:] = data_lower_type[i,sort2idx_l[i,:]];
                        data_upper_dist[i,:] = data_upper_dist[i,sort2idx_u[i,:]];
                        data_lower_dist[i,:] = data_lower_dist[i,sort2idx_l[i,:]];
                
                    dist_diff = np.zeros([nAtomPerLayer,nAtomPerLayer], dtype = float); # rate of difference on distance
                    type_ok = np.zeros([nAtomPerLayer,nAtomPerLayer], dtype = bool); # true if the type items are matching between a upper-atom and a lower-atom
                    for idx_upper in xrange(nAtomPerLayer):
                        for idx_lower in xrange(nAtomPerLayer):
                             dist_diff[idx_upper,idx_lower] = max(abs(np.divide(data_upper_dist[idx_upper,:]-data_lower_dist[idx_lower,:], data_upper_dist[idx_upper,:]+data_lower_dist[idx_lower,:])));
                             type_ok[idx_upper,idx_lower]= all(data_upper_type[idx_upper,:]==data_lower_type[idx_lower,:]);
                    minDiffRate = np.min(np.min(dist_diff), minDiffRate)

                    g = networkx.to_networkx_graph(dist_diff<=dist_tol_rate & type_ok); # 
                    # find the maximal matching of graph g
                    if len(networkx.maximal_matching(g))< nAtomPerLayer: 
                        polar=True;
                        break;
        
            if not(polar) and not(firstRound): # first round NonPolar
                firstRound = True
                layer_thickness = 0
                minNP_z = surface_z[-1]
                perodicity_lower_bound = abs(init_z - minNP_z); # the perodicity should be larger than this
            else:
                if not(polar) and firstRound and (abs(minNP_z[0] - surface_z[-1])>perodicity_lower_bound): # the second round non-polar
                    np.append(minNP_z, surface_z[-1]);
                    layer_thickness = layer_thickness + 1;
                    z_thickness = minNP_z[0]-minNP_z[-1];
        #             vz(3) = z_thickness;
        #             isNP = ismember(atom_coor(:,3), minNP_z);
        #             minNPStruc= atom_coor(isNP,:);
        #             minNPStruc(:,3)= minNPStruc(:,3) - minNP_z(end);
        #             atom_type = atom_type(isNP);
        #             atom_cell=[num2cell(minNPStruc) atom_type];
                    
        #             # create a new txt file that contains the maximal non-polar structure
        #             isWin = ~isempty(strfind(computer, 'PCWIN'));
        #             [pathstr,name,ext] = fileparts(filepath);
        #             mkdir(pathstr,'minNonPolar5');
        #             if isWin
        #                 pathstr=strcat(pathstr, '\minNonPolar5');
        #             else
        #                 pathstr=strcat(pathstr, '/minNonPolar5');
        #             end
        #             oldpath=cd(pathstr);
        #             fileID=fopen([name '-minNonPolar5' ext], 'w+');
        #             formatSpec1='lattice_vector \t %f \t %f \t %f \n';
        #             fprintf(fileID,formatSpec1,vx);
        #             fprintf(fileID,formatSpec1,vy);
        #             fprintf(fileID,formatSpec1,vz);
        #             formatSpec2='atom \t %f \t %f \t %f \t %s \n';
                    
        #             nrows= size(atom_cell,1);
        #             for row = 1:nrows:
        #                 fprintf(fileID,formatSpec2,atom_cell{row, :});
        #             fclose(fileID);
        #             #movefile([name '-nonpolar' ext], pathstr);
        #             cd(oldpath);
        #            
        #            top_z=minNP_z(1);
                    return polar, z_thickness, minDiffRate
        
            if polar and firstRound:
                layer_thickness = layer_thickness + 1;
                np.append(minNP_z,surface_z[-1]);

            data_delete = (atom_coor[:, 2]==surface_z[-1]);
            data_delete = data_delete[~forDelete];
            forDelete= (forDelete | (atom_coor[:, 2]==surface_z[-1]));
            atom_data_dist = atom_data_dist[data_delete,:];
            atom_data_dist = atom_data_dist[:,data_delete];
            atom_data_type = atom_data_type[data_delete,:];
            atom_data_type = atom_data_type[:,data_delete];
            surface_n=surface_n-1;
            surface_z = np.delete(surface_z, -1);

        
        z_thickness=float("nan");
        layer_thickness = float("nan");
        #top_z = float("nan");
        #perodicity_lower_bound = float("nan");

        return polar, z_thickness, minDiffRate




#########################################
#					#
#					#
#          Start program 		#
#					#
#					#
#########################################

#start timer
tStart = time.time()
atomTyp = {}
# Read input
subCIF, useMillerList, maxMillerInd, MillerIndList, nL, \
checkPolarity, z_digit, dis_tol_rate = readInput(inputFile)

# create a list of Miller indices
if not useMillerList:
	MillerList = createMillerList(maxMillerInd)
else:
	MillerList = MillerIndList
	print MillerList

#Read CIF file
print
print "********************************************************"
print "Reading structure from .cif file"
idMat,transM,atoms,positions,atomTyp=ReadCIF(subCIF,atomTyp)

# Construt big bulk material that will be reused in all calculations
print "Construction big bulk structure... This might take time."
bigBulk = Surface(transM,positions,atoms,atomTyp,np.array((0,0,0)))
bigBulk.bulkNEW(40)
print "Bulk structure complete"
print "********************************************************"
print 

primFailed = [] # stuctures for which cound find primtive vectors
notExist = [] # stuctures for which given orientation does not exist
# Start the loop on Miller indices 
isPolar = []; # if the structure is polar or not
perodicity  = []; # if polar, report the perodicity in angstrom; otherwise report NaN
minDiffRate = []; # the minimal difference rate
polarFilename=subCIF.split(".")[0]+"-poalrity.txt"
file = open(polarFilename, 'w+')
file.write("Filename\t\tisPolar\t\tperodicity\t\tminDiffRate\n")
file.close()

for subMillerString in MillerList:
	#MATERIAL 1
	print
	print "*****************************************"
	print "Constructing bulk Substrate"
	print "Orientation: %s"%subMillerString
	#idMat,transM,atoms,positions,atomtyp=ReadCIF(subCIF)
	Miller = getMillerFromString(subMillerString)
	#nbulk = max(Miller)*2 # use the max Miller index to ensure non-empty surface

	Sub = Surface(transM,positions,atoms,atomTyp,Miller)
#	Sub.bulk(nbulk)
	Sub.positions = bigBulk.positions.copy()
	Sub.atoms = bigBulk.atoms.copy()
	Sub.positionsSmall = bigBulk.positionsSmall.copy()
	Sub.atomsSmall = bigBulk.atomsSmall.copy()

	Sub.construct()
	Sub.plane()
	if not Sub.exists: 
		print "!!! Orientation does not exists!!!"
		print "!!! Proceeding to next one !!!"
		print "*****************************************"
		notExist.append(subMillerString)
		continue # if given plane does not exists, continue to next one

	Sub.initpvecNEW2()
	if not Sub.exists: 
		print "!!! Failed to find primitive vectors !!!"
		print "!!! Proceeding to next orientation  !!!"
		print "*****************************************"
		primFailed.append(subMillerString)
		continue # if given plane does not exists, continue to next one

	Sub.primitivecell()

	vecsS = Superlattice(Sub.a,Sub.b,1) # multiplication of unit cells
	# vecsS has set of vecors that span on the lattice. They give the same area,
	# but the two vectors might be of different lengths. Find such a pair, for
	# which length is similar, so we have the most "square" lattice
	maxRatio = 0
	ii = 0 #counter
	for vec in vecsS:
		normA = np.linalg.norm(vec[0])
		normB = np.linalg.norm(vec[1])
		ratio = normA/normB
		if normA >= normB : ratio = normB/normA
		if ratio > maxRatio:
			maxRatio = ratio
			vecIndex = ii
		ii += 1

	#Construct big planes for Substrate and Deposit
	#print "CREATING BULK - this may take time"
	#nbulk2 = max(Miller)*2
	#Sub.bulk(nbulk2)
	Sub.positions = bigBulk.positions.copy()
	Sub.atoms = bigBulk.atoms
	Sub.construct()
	Sub.plane()

	print "CREATING STRUCTURE"
	iface = Interface(vecsS[vecIndex],Sub,nL)

	print "OUTPUTTING COORDINATES"
	# Output in cartesian coordiantes
	millerLab= subMillerString.split()
	strufilename=subCIF.split(".")[0]+"-%s%s%s"%(millerLab[0],
			                            millerLab[1],
						    millerLab[2],)

	"""
	file = open(strufilename+".xyz",'w')
	natoms = len(iface.SSurfPos)
	file.write("%i\n\n"%natoms)
	ii = 0
	for i in iface.SSurfPos:
		file.write("%5s  %12.6f  %12.6f  %12.6f\n"%(iface.SSurfAtm[ii], i[0],i[1],i[2]))
		ii+= 1
	file.close()
	"""

	# Output in FHI-AIMS .in format
	file = open(strufilename+".in",'w')
	file.write("lattice_vector   %12.6f   %12.6f   %12.6f\n"%(iface.vecSubR[0][0],iface.vecSubR[0][1],iface.vecSubR[0][2]))
	file.write("lattice_vector   %12.6f   %12.6f   %12.6f\n"%(iface.vecSubR[1][0],iface.vecSubR[1][1],iface.vecSubR[1][2]))
	file.write("lattice_vector   %12.6f   %12.6f   %12.6f\n"%(0.0, 0.0, nL))
	ii = 0
	atom_coor = np.zeros([len(iface.SSurfPos), 3])
    atomLabels = np.zeros([len(iface.SSurfPos), 1])
    for i in iface.SSurfPos:
        file.write("atom %12.6f  %12.6f  %12.6f  %5s\n"%(i[0],i[1],i[2],atomTyp[iface.SSurfAtm[ii]]))
        atom_coor[ii,:] = i[range(3)];
        atomLabels[ii] = iface.SSurfAtm[ii];
        ii+= 1
    file.close()
    # check the polarity here
    if checkPolarity:
        print "Check the poalrity"
        polar, ped, minDiff = iface.checkPolar(atom_coor,atomLabels, z_digit, dis_tol_rate);
        isPolar.append(polar);
        perodicity.append(ped);
        minDiffRate.append(minDiff);
        if polar:
            print "The structure is polar with perodicity %s angstrom."%ped  
        else:
            print "The structure is non-polar, with a minimal distance difference rate %s."%minDiff 
        file = open(polarFilename, 'a')
        file.write(strufilename+"\t\t%d\t\t%f\t\t%f\n"%(polar,ped,minDiff))
        file.close()

	# end of the loop on Miller indices 

# nlayers? #
	print "*****************************************"

# Stop timer
tStop = time.time()

runtime = calcRuntime(tStart,tStop)

#Output statistics
nStruc = len(MillerList)
nPrimFailed = len(primFailed)
nNotExist = len(notExist)
nCreated = nStruc - nPrimFailed - nNotExist
file = open('stats.out','w')
file.write("Max. Miller Index: %i\n"%maxMillerInd)
file.write("Total number of possible orientations: %i\n"%nStruc)
file.write("Number created orientations: %i\n"%nCreated)
file.write("Number of not existing orientations: %i\n"%nNotExist)
for i in notExist:
	file.write("%s\n"%i)

file.write("\nNumber of orientations where primitive vectors failed: %i\n"%nPrimFailed)
for i in primFailed:
	file.write("%s\n"%i)
file.write("\nRuntime: %i min %i sec\n"%(runtime[0],runtime[1]))
file.close()

# End of program