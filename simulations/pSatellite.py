# Plasmid Inheritance Simulator V3
# Stratton Georgoulis
# Barrick lab
# Created 11/8/17
# Last edited 12/118/17
# In this version, a parent cell is chosen from the population with 20 plasmids. Plasmids replicate to 40 inside parent
# cell, then randomly segregatie into daughter cells so each daughter gets 20 plasmids. One daughter replaces the parent,
# and the other duaghter replaces a cell chosen at random. If one daughter only gets miniplasmid, it is possible for one of
# the two daughters to have a fitness of zero. If this is the case, the duaghter with non-zero fitness replaces the parent,
# but the zero-fitness daughter doesn't enter the population and no random cell is chosen. Miniplasmids also have a 
# slightly lower cost to the cell than previous versions. Cost of all plasmids reduced to half since cells now carry a 
# maximum of 20 to prevent cells having a fitness of lower than 0. This model deosn't include possible integration
# of the antibiotic resistance gene onto the chromosome.
import random
import time


# Cell population stored as dictionary of lists:
#   Keys are tuples as number of (FP, MP, DP) in that category of cell
#   Values are lists with (Fitness, Selection Probability, Cell Number)

#Global settings
plasmidsPerCell = 20

# Relativ fitnesses of plasmids for intra-cell replication
fullPlasmidFitness = 1   # Fitness for FP, MP, and DP. FP and DP arbitrarily set to 0.4
miniPlasmidFitness = 1.2 # 1.2 to 1.5 gives correct ratio
deletionPlasmidFitness = 1.05

FPtoMPrate = (1 * 10 ** -5) # Mutation rates for full plasmid to mini plasmid
FPtoDPrate = (1 * 10 ** -5) # Mutation rates for full plasmid to deletion plasmid

#Set to empty string to not create file
populationFileName = "population.csv"
summaryFileName = "summary.csv"

# Models for plasmid replication and segregation
# 1 = First, replicate plasmidsPerCell split new plasmids. Then split exactly 50/50 into daughter cells.
# 2 = First replicate plasmidsPerCell split new plasmids. Then, coin flip for which plasmids end up in each daughter cell.
# 4 = First, segregate plasmid via coin flip to new cells. Then, replicate each cell up to plasmidsPerCell split.
plasmidReplicationSegregationModel = 4


#Relative fitness
# all full plasmids     = 0.45
# contains mini plasmid = 0.81
# fraction mini plasmid = 0.8
# all deletion plasmid  = 0.90

fitnessAllFP = 0.45
fitnessContainsMP = 0.81
averageFractionMP = 0.8
fitnessAllDP = 0.9

fitnessCostFP = (1 - fitnessAllFP)/plasmidsPerCell
fitnessCostMP = ((1 - fitnessCostFP * (plasmidsPerCell * (1-averageFractionMP))) - fitnessContainsMP) / (plasmidsPerCell * averageFractionMP)
fitnessCostDP = (1 - fitnessAllDP)/plasmidsPerCell

print "\nFitness Model"
print "  Full Plasmid Cost = " + str(fitnessCostFP)
print "  Mini Plasmid Cost = " + str(fitnessCostMP)
print "  Del  Plasmid Cost = " + str(fitnessCostDP)
print ""


def computeFitness(cellType):
	if cellType[0] == 0 and cellType[2] == 0: # If the cell doesn't have antibiotic resistance gene
		fitness = 0 # Fitness is 0
	else:
		fitness = 1 # If an antibiotic resistance gene is present, fitness depends on plasmids in cell
		fitness -= fitnessCostFP * cellType[0] # Fitness cost for full plasmid
		fitness -= fitnessCostMP # Fitness cost for mini plasmid.
		fitness -= fitnessCostDP # Fitness cost for deletion plasmid
	return fitness # Return fitness of a cell containing those plasmids

def computeSelectionProb(states, pop):
	popFit = 0
	for key in states:
		states[key][1] = states[key][0] * states[key][2]
		popFit += states[key][1]
	for key in states:
		states[key][1] /= popFit
	return states

def pickCellToDivide(states):
	randNum = random.uniform(0, 1)
	s = 0
	for key in states:
		s += float(states[key][1])
		if s >= randNum:
			return key

def pickRandCellToReplace(states, pop):
	randNum = random.uniform(0, 1)
	s = 0
	#total = 0
	#for key in states:
	#	total += states[key][2]
	#replaced = False
	for key in states:
		s += float(states[key][2]) / pop
		#print str(key) + " " + str(s) + "\n"
		if s >= randNum:
			states[key][2] -= 1
			if (states[key][2] == 0):
				del states[key]
	#		replaced = True
			break
			
	#if not replaced:
	#	print "DID NOT REPLACE CELL" + str(pop) + " " + str(total) + " " + str(randNum) + " " + str(s) + "\n"

# cell is a list with three values
def replicatePlasmidsInCell(cell,numNewPlasmids):

	pFitTotal = 0 # Sum total fitness of plasmids in cells
	pFitTotal += cell[0] * fullPlasmidFitness 
	pFitTotal += cell[1] * miniPlasmidFitness
	pFitTotal += cell[2] * deletionPlasmidFitness

	if (pFitTotal == 0):
		return

	#print str(cell) + " " + str(pFitTotal) + "\n"

	FPfit = (float(fullPlasmidFitness) * cell[0]) / pFitTotal # Determine relative fitnesses of plasmids in cell
	MPfit = (float(miniPlasmidFitness) * cell[1]) / pFitTotal
	DPfit = (float(deletionPlasmidFitness) * cell[2]) / pFitTotal
	
	#print str(FPfit) + " " + str(MPfit) + " " + str(DPfit) + "\n"

	pTot = sum(cell[0:3])
	for plasmid in range(0,numNewPlasmids):
		randNum = random.uniform(0, 1)
		if DPfit > randNum:
			cell[2] += 1
		elif DPfit + MPfit > randNum:
			cell[1] += 1
		else:
			randNum = random.uniform(0, 1)
			if FPtoMPrate > randNum:
				cell[2] += 1
			elif FPtoMPrate + FPtoDPrate > randNum:
				cell[1] += 1
			else:
				cell[0] += 1		

def divide(cell, states, pop):

	#print "Dividing cell  : " + str(cell) + "\n"
	#print "         state : " + str(states[cell]) + "\n"
	#total = 0
	#for key in states:
	#	total += states[key][2]
	#print " Total cells before : " + str(total) + "\n"

	
	states[cell][2] -= 1 # Subtract the cell about to reproduce from the population
	if (states[cell][2] == 0):
		del states[cell]

	cell = list(cell) # Convert the tuple to a list

	#For methods that replicate BEFORE segregate
	if (plasmidReplicationSegregationModel == 1) or (plasmidReplicationSegregationModel == 2):
		replicatePlasmidsInCell(cell,plasmidsPerCell)

	daughter1 = [0, 0, 0] # Make two empty daughter cells
	daughter2 = [0, 0, 0]


### Divide plasmids
	plasmidList = []
	plasmidList += cell[0] * [0]
	plasmidList += cell[1] * [1]
	plasmidList += cell[2] * [2]
		
	# Equal
	if (plasmidReplicationSegregationModel == 1):
		random.shuffle(plasmidList)
		for p in range(0,plasmidsPerCell):
			daughter1[plasmidList[p]] += 1
		for p in range(plasmidsPerCell,2*plasmidsPerCell):
			daughter2[plasmidList[p]] += 1	
	# Coin flip
	elif (plasmidReplicationSegregationModel == 2) or (plasmidReplicationSegregationModel == 4):
		for p in plasmidList:
			randNum = random.uniform(0, 1)
			if randNum < 0.5:
				daughter1[p] += 1
			else:
				daughter2[p] += 1
	else:
		exit()
	
		#For methods that replicate AFTER segregate
	if (plasmidReplicationSegregationModel == 4):
		replicatePlasmidsInCell(daughter1,plasmidsPerCell-daughter1[0]-daughter1[1]-daughter1[2])
		replicatePlasmidsInCell(daughter2,plasmidsPerCell-daughter2[0]-daughter2[1]-daughter2[2])
	
	daughter1 = tuple(daughter1)
	daughter2 = tuple(daughter2)
	daughter1fit = computeFitness(daughter1)
	daughter2fit = computeFitness(daughter2)
	replaceRandCell = True # Switch to false if one of the daughters has fitness of 0

	if daughter1fit != 0:
		if daughter1 in states:
			states[daughter1][2] += 1 # If that cell type is already in states, add one to the count
		else:
			states[daughter1] = [daughter1fit, None, 1] # Otherwise, new key and compute fitness
#	else:
#		replaceRandCell = False # If fitness of daughter cell is 0, don't put in population or replace random cell


	if daughter2fit != 0:
		if daughter2 in states:
			states[daughter2][2] += 1 # If that cell type is already in states, add one to the count
		else:
			states[daughter2] = [daughter2fit, None, 1] # Otherwise, new key and compute fitness
#	else:
#		replaceRandCell = False 

	if replaceRandCell: # This will be True if both cells added to population, False otherwise.
		pickRandCellToReplace(states, pop+1)

	total = 0
	for key in states:
		total += states[key][2]
	#print " Total cells after : " + str(total) + "\n"

def main():
	random.seed(7)
	
	t = int( time.time() * 1000.0 )
	random.seed( ((t & 0xff000000) >> 24) +
		((t & 0x00ff0000) >>  8) +
		((t & 0x0000ff00) <<  8) +
		((t & 0x000000ff) << 24)   )
	
	pop = int(input("Enter population size: ")) # The user enters the number of cells to be simulated.
	states = {} # Create a dictionary of "cell states" that keeps track of cells in population
	summaryfile = open ("summary.txt", "w") # Open file to record results

	popfile = 0
	if populationFileName:
		popfile = open (populationFileName, "w")
		popfile.write("gen,FP,MP,DP,num,fitness\n")
		
	# In state dictionary:
	# Each key will be plasmid count: (full plasmids, mini plasmids, deletion plasmids, chromosomal integration)
	# Each value will be a list: [the computed fitness of cells with that plasmid count, fitness relative to population, and the number of cells of that type].

	states[(plasmidsPerCell, 0, 0)] = [computeFitness((plasmidsPerCell, 0, 0)), None, pop] # Create state dictionary, all cells have only full plasmids
	keepGoing = True # User will switch this to false to exit program
	generation = 0 # A counter for keeping track of number of generation simulated

	while keepGoing: 
		generations = int(input("Enter the number of generations to be simulated: ")) # User inputs number of generations to be simulated
		for gen in range(generations):
			generation += 1 # Add one to counter of number of generations simulated
			for eachCell in range(pop): # A number of cells equal to the size of the population will be chosen to divide each generation
				states = computeSelectionProb(states, pop)
				cellToDiv = pickCellToDivide(states)
				divide(cellToDiv, states, pop)

			FPonlyCount = 0 # Reset cell counters
			FPcount = 0
			MPcount = 0
			MPmajority = 0
			DPcount = 0
			DPonlyCount = 0
			INcount = 0
			INonlyCount = 0
			FPtotal = 0
			MPtotal = 0
			DPtotal = 0
			total = 0
			totalPlasmids = 0
			# Now, add up number of cells in categories above from dictionary. Key[2] for each entry
			# contains the number of cells of that type.

			for key in states:
				total += states[key][2]
				totalPlasmids += states[key][2] * (key[0] + key[1] + key[2]) # Total number of plasmids in the population at any given time

				if key[0] > 0 and key[1] == 0 and key[2] == 0: # Cells with only full plasmids
					FPonlyCount += states[key][2]
				if key[1] > 0: # Cells with at least one miniplasmid
					MPcount += states[key][2]
				if key[1] > (key[0] + key[1] + key[2]) * 0.5: # Cells with majority (6 or more) miniplasmdis
					MPmajority += 1
				if key[2] > 0: # Cells with at least one deletion plasmid
					DPcount += states[key][2]
				if key[2] > 0 and key[0] == 0 and key[1] == 0: # Cells with deletion plasmid but no full or mini plasmid
					DPonlyCount += states[key][2]

				FPtotal += key[0] * states[key][2] # Total number of full plasmids
				MPtotal += key[1] * states[key][2] # Total number of mini plasmids
				DPtotal += key[2] * states[key][2] # Total number of deletion plasmids
			FPratio = (float(FPtotal) / totalPlasmids) * 100 # Full plasmids as percentage of total plasmids
			MPratio = (float(MPtotal) / totalPlasmids) * 100 # Mini plasmids as percentage of total plasmids
			DPratio = (float(DPtotal) / totalPlasmids) * 100 # Deletion plasmids as percentage of total plasmds

			print ("Generation " + str(generation) + ":\n")
			print ("Total cells in population: " + str(total))
			print ("Number of cells with only full plasmids: " + str(FPonlyCount))
			print ("Number of cells containing at least one mini plasmid: " + str(MPcount))
			print ("Number of cells containing majority mini plasmids: " + str(MPmajority))
			print ("Number of cells containing at least one deletion plasmid: " + str(DPcount))
			print ("Number of cells containing only deletion plasmid: " + str(DPonlyCount))
			print ("Total plasmids in population: " + str(totalPlasmids))
			print ("Full plasmids as percent of total Plasmids: " + str(FPratio) + "%")
			print ("Mini plasmids as percent of total plasmids: " + str(MPratio) + "%")
			print ("Deletion plasmids as percent of total plasmids: " + str(DPratio) + "%")
			print ("__________________________________________________________________________\n")

			# Also write the results to a .txt so they can be referenced later

			summaryfile.write (("Generation " + str(generation) + ":"))
			summaryfile.write (("Number of cells with only full plasmids: " + str(FPonlyCount) + "\n"))
			summaryfile.write (("Number of cells containing at least one mini plasmid: " + str(MPcount) + "\n"))
			summaryfile.write (("Number of cells containing majority mini plasmids: " + str(MPmajority) + "\n"))
			summaryfile.write (("Number of cells containing at least one deletion plasmid: " + str(DPcount) + "\n"))
			summaryfile.write (("Number of cells containing only deletion plasmid: " + str(DPonlyCount) + "\n"))
			summaryfile.write (("Full plasmids as percent of total Plasmids: " + str(FPratio) + "%\n"))
			summaryfile.write (("Mini plasmids as percent of total plasmids: " + str(MPratio) + "%\n"))
			summaryfile.write (("Deletion plasmids as percent of total plasmids: " + str(DPratio) + "%\n"))
			summaryfile.write ("__________________________________________________________________________\n\n")

			if populationFileName:
				for key in states:
					popfile.write(str(generation) + "," + str(key) + "," + str(states[key][2]) + "," + str(states[key][0]) + "\n")
			
		command = input("Would you like to continue the simulation? y/n: ") # Ask user whether to continue simulation
		while command != "y" and command != "n":
			command = input("Would you like to continue? y/n: ") # If invalid input, prompt the user again
		if command == "y": # Continue the simulation
			keepGoing = True 
		else:
			keepGoing = False # Breaks loop to exit simulation
	outfile.close()
	print("Exiting...")
main()
