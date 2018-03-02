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
print("Would you like to run simulation or run controls?")
runControls = input("Type s for simulation or c for controls: ")
while runControls != "s" and runControls != "c":
	runControls = input("Please type a single character, s or c: ")
if runControls == "c":
	runControls = True
else:
	runControls = False

if runControls == True:
	print ("Enter desired control to run: ")
	print ("For controls 1-4, all mutation rates are set to 0: ")
	print ("1: 90% of starting population has only full plasmids, while 10% has only deletion plasmids.")
	print ("2: 50% of starting population has only full plasmids, and 50% has only satellite plasmids.")
	print ("3: 90% of starting population has only full plasmids, 10% has AbR integration in chromosome and no plasmids.")
	print ("4: 90% of starting population has only full plasmids, and 10% has AbR integration in chromosome and a full set of full plasmids.")
	print ("5: Population starts with only full plasmid, FP --> MP and INT rates are 0.")
	initialPopulation = int(input("Enter desired control number: "))
else:
	initialPopulation = 0 # Indicates that no control is being run

print("Would you like to use default values for mutation rates, fitnesses, etc.?")
useDefaults = input("y/n: ")
while useDefaults != "y" and useDefaults != "n":
	useDefaults = input("Please type a single character, y/n: ")
if useDefaults == "y":
	useDefaults = True
else:
	useDefaults = False

if useDefaults:
	plasmidsPerCell = 20
	# Relative fitnesses of plasmids for intra-cell replication
	fullPlasmidFitness = 1
	miniPlasmidFitness = 1.2
	deletionPlasmidFitness = 1.05
	# Mutation rates for full plasmid to convert to miniplasmid or deletion plasmid
	FPtoMPrate = 3 * (10 ** -5)
	FPtoDPrate = 10 ** -5
	IntRate = 10 ** -5
	# Relative fitness
	# all full plasmids     = 0.45
	# contains mini plasmid = 0.81
	# fraction mini plasmid = 0.8
	# all deletion plasmid  = 0.90
	fitnessAllFP = 0.45
	fitnessContainsMP = 0.81
	averageFractionMP = 0.80
	fitnessAllDP = 0.90
	fitnessInt = 0.90
else:
	plasmidsPerCell = int(input("Enter number of plasmids per cell: "))
	# Relative fitnesses of plasmids for intra-cell replication
	fullPlasmidFitness = float(input("Enter full plasmid fitness: "))
	miniPlasmidFitness = float(input("Enter miniplasmid fitness: "))
	deletionPlasmidFitness = float(input("Enter deletion plasmid fitness: "))
	fitnessInt = float(input("Enter fitness for AbR integration: "))
	# Mutation rates for full plasmid to convert to miniplasmid or deletion plasmid
	FPtoMPrate = float(eval(input("Enter mutation rate for full plasmid to miniplasmid: ")))
	FPtoDPrate = float(eval(input("Enter mutation rate for full plasmid to deletion plasmid: ")))
	IntRate = float(eval(input("Enter rate of AbR integration to chromosome: ")))
	
	fitnessAllFP = float(input("Enter relative fitness for all full plasmids: "))
	fitnessContainsMP =  float(input("Enter relative fitness for cells with miniplasmid: "))
	averageFractionMP = float(input("Enter average fraction on miniplasmid: "))
	fitnessAllDP =  float(input("Enter relative fitness for all deletion plasmids: "))

if initialPopulation in (1, 2, 3, 4):
	print("ping")
	FPtoMPrate = 0
	FPtoDPrate = 0
	IntRate = 0
elif initialPopulation == 5:
	FPtoMPrate = 0
	IntRate = 0

fitnessCostFP = (1 - fitnessAllFP)/plasmidsPerCell
fitnessCostMP = ((1 - fitnessCostFP * (plasmidsPerCell * (1-averageFractionMP))) - fitnessContainsMP) / (plasmidsPerCell * averageFractionMP)
fitnessCostDP = (1 - fitnessAllDP)/plasmidsPerCell
fitnessCostInt = 1 - fitnessInt

print ("\nFitness Model")
print ("  Full Plasmid Cost = " + str(fitnessCostFP))
print ("  Mini Plasmid Cost = " + str(fitnessCostMP))
print ("  Del  Plasmid Cost = " + str(fitnessCostDP))
print ("  Integration Cost = " + str(fitnessCostInt))

#Set to empty string to not create file
populationFileName = "population.csv"
summaryFileName = "summary.csv"

print("\nModels for plasmid replication and segregation. In all models, the parent first has a chance of chromosomal integration of AbR; if AbR present, always passed to daughter cells.")
print("1 = First, replicate plasmidsPerCell split new plasmids. Then split exactly 50/50 into daughter cells.")
print("2 = First replicate plasmidsPerCell split new plasmids. Then, coin flip for which plasmids end up in each daughter cell.")
print("4 = First, segregate plasmid via coin flip to new cells. Then, replicate each cell up to plasmidsPerCell split.\n")
plasmidReplicationSegregationModel = int(input("Enter desired plasmid segregation model: "))

print("Daughter cells can either be replaced in initial population, or can be put in new population each generation.")
print("1 = Daughter cells are replaced in initial population. Number of divisions in a generation equals the size of the population.")
print("2 = Daughter cells placed in offspring population, which becomes next generation after reaching number of cells in ancestral population.")
daughterPlacement = int(input("Enter number for desired model: "))


def computeFitness(cellType):
	if cellType[0] == 0 and cellType[2] == 0 and cellType[3] == 0: # If the cell doesn't have antibiotic resistance gene
		fitness = 0 # Fitness is 0
	else:
		fitness = 1 # If an antibiotic resistance gene is present, fitness depends on plasmids in cell
		fitness -= fitnessCostFP * cellType[0] # Fitness cost for full plasmid
		fitness -= fitnessCostMP * cellType[1] # Fitness cost for mini plasmid.
		fitness -= fitnessCostDP * cellType[2]# Fitness cost for deletion plasmid
		fitness -= fitnessCostInt * cellType[3] # Fitness cost for AbR integration in chromosome
	return fitness # Return fitness of a cell containing those plasmids

def computeSelectionProb(states):
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

	pTot = cell[0] + cell[1] + cell[2]
	for plasmid in range(0,numNewPlasmids):
		randNum = random.uniform(0, 1)
		if DPfit > randNum:
			cell[2] += 1
		elif DPfit + MPfit > randNum:
			cell[1] += 1
		else:
			randNum = random.uniform(0, 1)
			if FPtoMPrate > randNum:
				cell[1] += 1
			elif FPtoMPrate + FPtoDPrate > randNum:
				cell[2] += 1
			else:
				cell[0] += 1	

def initializePopulation (pop, states, initialPopulation):
	if initialPopulation == 1:
		states[(plasmidsPerCell, 0, 0, 0)] = [computeFitness((plasmidsPerCell, 0, 0, 0)), None, int(0.9 * pop)]
		states[(0, 0, plasmidsPerCell, 0)] = [computeFitness((0, 0, plasmidsPerCell, 0)), None, int(0.1 * pop)]
	elif initialPopulation == 2:
		states[(plasmidsPerCell, 0, 0, 0)] = [computeFitness((plasmidsPerCell, 0, 0, 0)), None, int(0.5 * pop)]
		states[(0, plasmidsPerCell, 0, 0)] = [computeFitness((plasmidsPerCell, 0, 0, 0)), None, int(0.5 * pop)]
	elif initialPopulation == 3:
		states[(plasmidsPerCell, 0, 0, 0)] = [computeFitness((plasmidsPerCell, 0, 0, 0)), None, int(0.9 * pop)]
		states[(0, 0, 0, 1)] = [computeFitness((plasmidsPerCell, 0, 0, 0)), None, int(0.1 * pop)]
	elif initialPopulation == 4:
		states[(plasmidsPerCell, 0, 0, 0)] = [computeFitness((plasmidsPerCell, 0, 0, 0)), None, int(0.9 * pop)]
		states[(plasmidsPerCell, 0, 0, 1)] = [computeFitness((plasmidsPerCell, 0, 0, 0)), None, int(0.1 * pop)]
	else:
		states[(plasmidsPerCell, 0, 0, 0)] = [computeFitness((plasmidsPerCell, 0, 0, 0)), None, pop] # all cells have only full plasmids

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

	randNum = random.uniform(0, 1)
	if IntRate > randNum:
		cell[3] = 1

	daughter1 = [0, 0, 0, 0] # Make two empty daughter cells
	daughter2 = [0, 0, 0, 0]

	if cell[3] == 1: # If the parent has AbR integration, so will both daughters.
		daughter1 = [0, 0, 0, 1]
		daughter2 = [0, 0, 0, 1]

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

def divideIntoNextGeneration(cell, states, nextGeneration):
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

	plasmidsCanReplicate = True
	if cell[0] == 0 and cell[2] == 0:  # If cell lacks full and deletion plasmids, plasmids can't replicate
		plasmidsCanReplicate = False

	#For methods that replicate BEFORE segregate
	if ((plasmidReplicationSegregationModel == 1) or (plasmidReplicationSegregationModel == 2)) and plasmidsCanReplicate:
		replicatePlasmidsInCell(cell,plasmidsPerCell)

	randNum = random.uniform(0, 1)
	if IntRate > randNum:
		cell[3] = 1

	daughter1 = [0, 0, 0, 0] # Make two empty daughter cells
	daughter2 = [0, 0, 0, 0]

	if cell[3] == 1:
		daughter1 = [0, 0, 0, 1]
		daughter2 = [0, 0, 0, 1]

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
	if (plasmidReplicationSegregationModel == 4) and plasmidsCanReplicate:
		replicatePlasmidsInCell(daughter1,plasmidsPerCell-daughter1[0]-daughter1[1]-daughter1[2])
		replicatePlasmidsInCell(daughter2,plasmidsPerCell-daughter2[0]-daughter2[1]-daughter2[2])
	
	daughter1 = tuple(daughter1)
	daughter2 = tuple(daughter2)
	daughters = [daughter1, daughter2]

	for daughter in daughters:
		if daughter in nextGeneration:
			nextGeneration[daughter][2] += 1
		else:
			nextGeneration[daughter] = [computeFitness(daughter), None, 1]

	# total = 0
	# for key in states:
	# 	total += states[key][2]
	#print " Total cells after : " + str(total) + "\n"



def main():
	seednum = int(input("Enter an integer value for seed: "))
	random.seed(seednum)
	
	# t = int( time.time() * 1000.0 )
	# random.seed( ((t & 0xff000000) >> 24) +
	# 	((t & 0x00ff0000) >>  8) +
	# 	((t & 0x0000ff00) <<  8) +
	# 	((t & 0x000000ff) << 24)   )
	
	pop = int(input("Enter population size: ")) # The user enters the number of cells per generation to be simulated.
	divisionsPerGeneration = pop
	if daughterPlacement == 2:
		divisionsPerGeneration /= 2
		divisionsPerGeneration = int(divisionsPerGeneration)
	states = {} # Create a dictionary of "cell states" that keeps track of cells in population
	summaryfile = open ("summary.txt", "w") # Open file to record results

	popfile = 0
	if populationFileName:
		popfile = open (populationFileName, "w")
		popfile.write("gen,FP,MP,DP,num,fitness\n")
		
	# In state dictionary:
	# Each key will be plasmid count: (full plasmids, mini plasmids, deletion plasmids, chromosomal integration)
	# Each value will be a list: [the computed fitness of cells with that plasmid count, fitness relative to population, and the number of cells of that type].

	print("initialPopulation: ", initialPopulation)

	initializePopulation (pop, states, initialPopulation)
	keepGoing = True # User will switch this to false to exit program
	generation = 0 # A counter for keeping track of number of generation simulated

	while keepGoing: 
		generations = int(input("Enter the number of generations to be simulated: ")) # User inputs number of generations to be simulated
		for gen in range(generations):
			generation += 1 # Add one to counter of number of generations simulated
			if daughterPlacement == 2:
				nextGeneration = {}
			for eachCell in range(divisionsPerGeneration): # A number of cells equal to the size of the population will be chosen to divide each generation
				states = computeSelectionProb(states)
				cellToDiv = pickCellToDivide(states)
				if daughterPlacement == 1:
					divide(cellToDiv, states, pop)
				elif daughterPlacement == 2:
					divideIntoNextGeneration(cellToDiv, states, nextGeneration)
			if daughterPlacement == 2:
				states = nextGeneration
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
				if key[1] > (key[0] + key[1] + key[2]) * 0.5: # Cells with majority miniplasmdis
					MPmajority += states[key][2]
				if key[2] > 0: # Cells with at least one deletion plasmid
					DPcount += states[key][2]
				if key[2] > 0 and key[0] == 0 and key[1] == 0: # Cells with deletion plasmid but no full or mini plasmid
					DPonlyCount += states[key][2]
				if key[3] == 1:  # Cells with AbR integration
					INcount += states[key][2]
				if key[0] == 0 and key[1] == 0 and key[2] == 0 and key[3] == 0: # Cells with AbR integration and no other plasmids
					INonlyCount += states[key][2]

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
			print ("Number of cells with AbR intregrated into chromosome: " + str(INcount))
			print ("Number of cells with AbR integraged into chromosome with no plasmids: " + str(INonlyCount))
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
			summaryfile.write (("Number of cells with AbR intregrated into chromosome: " + str(INcount)))
			summaryfile.write (("Number of cells with AbR integraged into chromosome with no plasmids: " + str(INonlyCount)))
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
	summaryfile.close()
	print("Exiting...")
main()