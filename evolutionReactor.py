
import sys, math, random, os, matplotlib.pyplot as plt
import numpy as np

class Cell:
    def __init__(self, seq_1, struc_1 ):
        self.sequence_1 = seq_1
        self.structure_1 = struc_1

# compute structure using RNAfold
def computeStruc( sequence):
    cmd = "echo "+ sequence + " > tmpq6; RNAfold < tmpq6 | tail -n1 | cut -d' ' -f1 > tmpq6out"
    #res = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    os.system( cmd )
        
    for line in open( 'tmpq6out', 'r'):
        res = line
        res = res.replace('\n', '')
        break
        
    return res

# compute distance using RNAdistance
def bp_distance( fold, target ):
    fout = open('tmp2q6', 'w')
    fout.write(fold + '\n' + target)
    fout.close()
        
    cmd = "RNAdistance < tmp2q6 | cut -d' ' -f2  > tmp2q6out"
    #print (cmd)
    #sys.exit()
    os.system( cmd )
        
    res = ''
        
    for line in open( 'tmp2q6out', 'r'):
        res = line
        res = res.replace('\n', '')
        break
        
    #print ( fold, target, res )
    return int(res)


def populate(target, pop_size):
    
    population = []
    
    for i in range(pop_size):
        #get a random sequence to start with
        sequence = "".join([random.choice("AUCG") for _ in range(len(target))])
        #use nussinov to get the secondary structure for the sequence
        structure = computeStruc(sequence)
        #add a new Cell object to the population list, both chromosomes will be identical to start
        new_cell = Cell(sequence, structure)
        new_cell.id = i
        new_cell.parent = i
        population.append(new_cell)
    
    return population

def compute_fitness(population, target, beta):
    """
        Assigns a fitness and bp_distance value to each cell in the population.
        """
    #store the fitness values of each cell
    tot = []
    #iterate through each cell
    for cell in population:
        
        #calculate the bp_distance of each chromosome using the cell's structure
        bp_distance_1 = bp_distance(cell.structure_1, target)
        #bp_distance_2 = bp_distance(cell.structure_2[-1], target)
        
        #use the bp_distances and the above fitness equation to calculate the fitness of each chromosome
        fitness_1 = math.exp(beta * float(bp_distance_1) / len(target))
        #fitness_2 = math.exp(beta * float(bp_distance_2) / len(target))
        
        #get the fitness of the whole cell by multiplying the fitnesses of each chromosome
        cell.fitness = fitness_1 #* fitness_2
        
        #store the bp_distance of each chromosome.
        cell.bp_distance_1 = bp_distance_1
        #cell.bp_distance_2 = bp_distance_2
        
        
        #add the cell's fitness value to the list of all fitness values (used for normalization later)
        tot.append(cell.fitness)
    
    #normalization factor is sum of all fitness values in population
    norm = np.sum(tot)
    #divide all fitness values by the normalization factor.
    for cell in population:
        cell.fitness = cell.fitness / norm
    
    return None


def mutate(sequence, mutation_rate):
    """Takes a sequence and mutates bases with probability mutation_rate"""
    
    #start an empty string to store the mutated sequence
    new_sequence = ""
    #boolean storing whether or not the sequence got mutated
    mutated = False
    #go through every bp in the sequence
    for bp in sequence:
        #generate a random number between 0 and 1
        r = random.random()
        #if r is below mutation rate, introduce a mutation
        if r < mutation_rate:
            #add a randomly sampled nucleotide to the new sequence
            new_sequence = new_sequence + random.choice("AUCG")
            mutated = True
        else:
            #if the mutation condition did not get met, copy the current bp to the new sequence
            #print(bp)
            new_sequence = new_sequence + bp

    return (new_sequence, mutated)


def selection(population, target, mutation_rate, beta):
    """
        Returns a new population with offspring of the input population
        """
    
    #select the sequences that will be 'parents' and contribute to the next generation
    #look at the documentation for np.random.choice and its optional argument p
    parents = np.random.choice(population, len(population), p=[cell.fitness for cell in population], replace=True)
    
    #build the next generation using the parents list
    next_generation = []
    for i, p in enumerate(parents):
        new_cell = Cell(p.sequence_1, p.structure_1)
        new_cell.id = i
        new_cell.parent = p.id
        
        next_generation.append(new_cell)
    
    #introduce mutations in next_generation sequeneces and re-fold when a mutation occurs
    for rna in next_generation:
        mutated_sequence_1, mutated_1 = mutate(rna.sequence_1, mutation_rate)
        #mutated_sequence_2, mutated_2 = mutate(rna.sequence_2)
        
        #if mutation occured assign and fold the new sequence
        if mutated_1:
            rna.sequence_1 = mutated_sequence_1
            rna.structure_1 = computeStruc(mutated_sequence_1)
        #if mutated_2:
            #rna.sequence_2 = mutated_sequence_2
            #rna.structure_2 = nussinov(mutated_sequence_2)
        else:
            continue

    #update fitness values for the new generation
    compute_fitness(next_generation, target, beta)
    
    return next_generation


'''
def record_stats(pop, population_stats):
    """
        Takes a population list and a dictionary and updates it with stats on the population.
        """
    #print ( pop, 'pop length', len(pop) )
    for rna in pop:
        #print ( rna.bp_distance_1 )
        generation_bp_distance_1 = rna.bp_distance_1 #for rna in pop]


    #generation_bp_distance_1 = [rna.bp_distance_1 for rna in pop]
    #generation_bp_distance_2 = [rna.bp_distance_2 for rna in pop]
    
    mean_bp_distance_1 = np.mean(generation_bp_distance_1)
    #mean_bp_distance_2 = np.mean(generation_bp_distance_2)
    
    mean_fitness = np.mean([rna.fitness for rna in pop])
    
    
    population_stats.setdefault('mean_bp_distance_1', []).append(mean_bp_distance_1)
    #population_stats.setdefault('mean_bp_distance_2', []).append(mean_bp_distance_2)
    
    #population_stats.setdefault('mean_fitness', []).append(mean_fitness)
    
    return None
'''

def evolve(target, generations, pop_size, mutation_rate, beta):
    """
        Takes target structure and sets up initial population, performs selection and iterates for desired generations.
        """
    #store list of all populations throughotu generations [[cells from generation 1], [cells from gen. 2]...]
    populations = []
    #start a dictionary that will hold some stats on the populations.
    population_stats = {}
    
    #get a starting population
    initial_population = populate(target, pop_size)
    #compute fitness of initial population
    compute_fitness(initial_population, target, beta)
    
    #set current_generation to initial population.
    current_generation = initial_population
    
    
    #for p in current_generation:
    #    print ( p.sequence_1 )
    
    #sys.exit()
    
    
    #iterate the selection process over the desired number of generations
    for i in range(generations):
        
        #let's get some stats on the structures in the populations
        #record_stats(current_generation, population_stats)
        #print (i, '\n', current_generation )
        bpd = []
        for p in current_generation:
            bpd.append( p.bp_distance_1 )

        population_stats.setdefault('mean_bp_distance_1', []).append( np.mean( bpd ) )


        #add the current generation to our list of populations.
        populations.append(current_generation)
        
        #select the next generation
        new_gen = selection(current_generation, target, mutation_rate, beta)
        #set current generation to be the generation we just obtained.
        current_generation = new_gen
    
    return (populations, population_stats)


targets = [ '((((((((....))))))))', '((((..(((....)))))))', '(((....)))(((....)))']

mus = [ 0.01, 0.02, 0.05, 0.1 ]

colors = ['red', 'blue', 'green', 'yellow' ]

generations = 500
pop_size = 100

beta = -2

for trange in range(0, len(targets)):
    

    target = targets[ trange ]

    plotname = 'target' + str(trange) + '.png'

    plt.figure()
    
    for murange in range( len(mus) ):
        
        mu = mus[ murange ]
        
        print ( 'target: ', trange, ' mu: ', mu )
        
        mutation_rate = mu

        pops, pops_stats = evolve( target , generations, pop_size, mutation_rate, beta)
        
        pops_stats = pops_stats['mean_bp_distance_1']
        #print ( len(pops_stats))

        plt.scatter( range(0,generations), pops_stats, color=colors[murange], label = 'mu = '+str(mu) )

        #break

    plt.legend( loc='upper right')

    plt.xlabel('Generations')
    plt.ylabel('Average Distance')

    plt.title( 'target = ' + target )

    plt.savefig( plotname )

    #break


#print ( pops_stats )




















