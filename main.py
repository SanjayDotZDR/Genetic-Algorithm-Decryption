'''
Sanjay Suthakaran
June 18th, 2025
'''
import numpy as np
import random, string
import evaluation as eval




class GeneticAlgorithm:

    #population array of strings
   
    #TODO implement another typa crossover, clean up code and prepare some way to collect data for report.

    def __init__(self, filename, pop_size, max_gens, elitism_size, crossover_rate, mutation_rate, crossover_type, print_gen, data_filename):

        self.key_length, self.cipher_text = read_file(filename)
        self.pop_size = pop_size
        self.max_gens = max_gens
        self.elitism_size = elitism_size
        self.crossover_rate = crossover_rate
        self.mutation_rate = mutation_rate
        self.crossover_type = crossover_type
        self.print_gen = print_gen
        self.data_filename = data_filename


        #debug
        '''self.pop_size = 500
        self.max_gens = 4000
        self.elitism_size = 4
        self.crossover_rate = 0.7
        self.mutation_rate = 0.2
        self.key_length = 40'''

        

        #initialize population, will recieve an array of pop size
        #format (chromosome, fitness)

        self.best_solution = self.genetic_algorithm()
        
        print(str(self.best_solution) + " is the key!")

        
        
        


        
        #initialize pop (call the below func)

        #we will now have a population of size (pop_size). 


        #while term cond not met (will decide these later)

            #evaluate fitness for each chromosome, make sure to check fitness score meets target
            #put em into an array corresponding to each solution


        #new pop extend - select k most elite from fitness scores (sort it) and bring em to next gen. 

        #create new population----
        #tourney selection


        
        #print(self.pop_array)

    def genetic_algorithm(self):

        sol_found = False
        self.best_solution = None
        current_gen = 0
        no_change_count = 0
        prev_gen = None

        self.pop_tuple = self.create_initial_population()




        print("\n------ INITIAL POP --------")
        print(self.pop_tuple)
        print("\n--------------")

        #for x generations.... [ MAIN LOOP ]-------------------------------------------------------------------------
        while (current_gen < self.max_gens):
            

            #print("\n-----CURRENT GEN " + str(current_gen))
            #print(self.pop_tuple)

            #print("\n-------testestestest")


            #check fitness of each individual in population and append it to the tuple
            for i in range(self.pop_size):
                
                self.chr = self.pop_tuple[i][0]

                fitness = eval.fitness(self.chr, self.cipher_text)

                self.pop_tuple[i] = (self.chr,float(fitness))




            #check if target fitness reached (0.05)

            for self.chr in self.pop_tuple:
                
                if (self.chr[1] <= 0.05):

                    self.best_solution = (self.chr[0],self.chr[1])

                    sol_found = True

                    break
            if sol_found:
                break


            


            #termination condition not met, move onto creating next gen.
            #grab k most elite for next generation
            self.elite_tuple = self.choose_elite()

            #use tournament selection to grab popsize - elite individuals for next generation
            self.tournament_tuple = self.tournament_selection()
            
            #tuples in selected_parents to reproduce
            self.selected_parents = self.elite_tuple + self.tournament_tuple

            #uniform order crossover here
            if (self.crossover_type == "a"):
                
                self.pop_tuple = self.uniform_order_crossover()

            #two point crossover here
            else:
                
                self.pop_tuple = self.two_point_crossover()
                

            current_gen += 1

            #print every x gens
            if current_gen % self.print_gen == 0:
                
                for i in range(self.pop_size):
                    self.chr = self.pop_tuple[i][0]
                    fitness = eval.fitness(self.chr, self.cipher_text)
                    self.pop_tuple[i] = (self.chr,float(fitness))

                self.pop_tuple.sort(key=lambda x: x[1])
                print("\n\nCurrent Generation: " + str(current_gen) + "\nBest key of this Generation: " + str(self.pop_tuple[0][0]) + "\nFitness: " + str(self.pop_tuple[0][1]))

                self.avg = str(self.average_fitness())

                self.write_to_file(current_gen, str(self.pop_tuple[0][1]))







        
        print("\n--------")
        print("final gen:" + str(current_gen))

        for i in range(self.pop_size):
                
            self.chr = self.pop_tuple[i][0]
            fitness = eval.fitness(self.chr, self.cipher_text)
            self.pop_tuple[i] = (self.chr,float(fitness))


        for item in self.pop_tuple:
                print(f"Chromosome: {item[0]}, Fitness: {item[1]:.4f}")

        
        self.pop_tuple.sort(key=lambda x: x[1])

        print("\nBest attempt: \nKey: " + str(self.pop_tuple[0][0]) + "\nFitness: " + str(self.pop_tuple[0][1]))

        print(print("\nDecrypted: \n" + str(eval.decrypt(str(self.pop_tuple[0]),self.cipher_text))))

        return str(self.pop_tuple[0][0])
    

    def create_initial_population(self):
        #pop size is number of chromosomes
        #key length is the length of each chromosome

        avaliable_gene_pool = string.ascii_lowercase + "-"

        #temp tuple list
        self.temp = []


        for _ in range(self.pop_size):

            #create a random chromosome by randomly picking from gene pool and concatenating to string K times (key length)
            chr = ''.join(random.choices(avaliable_gene_pool, k=self.key_length))

            #append new chromosome to population
            self.temp.append((chr,0))

        return self.temp


    def choose_elite(self):
        
        #sort tuple from lowest to highest fitness.
        self.pop_tuple.sort(key=lambda x: x[1])
        
        temp_tuple = []

        for i in range(self.elitism_size):
            temp_tuple.append((self.pop_tuple[i][0], self.pop_tuple[i][1]))
            self.pop_tuple[i]

        

        return temp_tuple

    
    def tournament_selection(self):

        temp_tuple = []
        #repeat until popsize - elite individuals picked.
        for _ in range(self.pop_size - self.elitism_size):
            
            
            tournament_tuple = []

            k = random.randint(2,5)
            

            #pick k random individuals
            for _ in range(k):
                tournament_tuple.append(random.choice(self.pop_tuple))
            
            
            #compare the tournament contenders and pick best fitness
            tournament_tuple.sort(key=lambda x: x[1])
            chosen = tournament_tuple[0]

            temp_tuple.append(chosen)

        return temp_tuple


    def uniform_order_crossover(self):

        #uses selected parents tuple.

        temp_tuple = []

        for i in range(len(self.selected_parents) - 1 ):

            #every second element
            if i % 2 == 0:
                
                parent1 = self.selected_parents[i][0]
                parent2 = self.selected_parents[i + 1][0]


                #do a crossover
                if random.random() < self.crossover_rate:

                    mask = [random.randint(0,1) for _ in range(len(parent1))]


                    child1 = ''
                    child2 = ''

                    for i in range(len(parent1)):

                        #grab from "main" parent
                        if mask[i] == 1:
                            child1 += parent1[i]
                            child2 += parent2[i]

                        else:
                            child1 += parent2[i]
                            child2 += parent1[i]

                    #mutate
                    child1 = self.mutate(child1)
                    child2 = self.mutate(child2)

                    temp_tuple.append((child1,0))
                    temp_tuple.append((child2,0))
                    

                #skip crossover for these parents
                else:
                    temp_tuple.append((parent1,0))
                    temp_tuple.append((parent2,0))

        return temp_tuple


    def two_point_crossover(self):

         #uses selected parents tuple.

        temp_tuple = []

        for i in range(len(self.selected_parents) - 1 ):

            #every second element
            if i % 2 == 0:
                
                parent1 = self.selected_parents[i][0]
                parent2 = self.selected_parents[i + 1][0]


                #do a crossover only if crossover probability works + length of parent greater than 2 (2 point will not work otherwise)
                if (random.random() < self.crossover_rate) and (len(parent1) > 2):

                    #grab points randomly
                    p1 = random.randint(0, len(parent1) - 2)
                    p2 = random.randint(p1 + 1, len(parent1) - 1)

                    #string together children, using string splicing
                    #i.e take everything in parent 1 before p1 + take everything in parent 2 between p1 and p2 + take everything in parent1 after p2 = child1
                    child1 = parent1[:p1] + parent2[p1:p2] + parent1[p2:]
                    child2 = parent2[:p1] + parent1[p1:p2] + parent2[p2:]

                    #mutate
                    child1 = self.mutate(child1)
                    child2 = self.mutate(child2)

                    temp_tuple.append((child1,0))
                    temp_tuple.append((child2,0))
                    

                #skip crossover for these parents
                else:
                    temp_tuple.append((parent1,0))
                    temp_tuple.append((parent2,0))

        return temp_tuple


    def average_fitness(self):
     
        #list of all fitness scores in population tuple
        fitness_scores = [fitness for _, fitness in self.pop_tuple]


        return sum(fitness_scores) / len(fitness_scores)


    def write_to_file(self, current_gen, current_best):


        #writes in format Gen,Average fitness,Best fitness
        with open(f"Assignment1\{self.data_filename}", "a") as file:
            file.write(str(current_gen) + "," + self.avg + "," + str(current_best) + "\n")


    

    '''
    Reciprocol Exchange
    '''
    def mutate(self, individual):
        
        if random.random() < self.mutation_rate:
            
            genes = list(individual)
            gene1, gene2 = random.sample(range(len(genes)),2)


            swap1 = genes[gene1]
            swap2 = genes[gene2]

            genes[gene1] = swap2
            genes[gene2] = swap1

            return ''.join(genes)
            

        else:
            return individual  

        

def create_next_generation(fitness_scores):
    pass        


    


        

    def evaluate_population():
        pass

    def run(max_generations):
        pass






'''
Selection and Reproduction
'''









'''
Utility Functions and Helper
'''
def read_encrypted_data(filename):
    pass

def run_ga_on_file(filename):
    pass






    pass
def calculate_diversity(population):
    pass

def print_progress(generation, best_fitness):
    pass

def save_results(best_key, filename):
    pass




def read_file(filename):

    with open(filename, 'r') as file:

        lines = file.readlines()

        key_length = int(lines[0].strip())

        #skips first line and reads rest of file
        cipher_text = ''.join(line.strip() for line in lines[1:])
        
    
        return key_length, cipher_text






if __name__ == "__main__":
    print("\n-Genetic Algorithm Assignment-")

    filename = input("\nWhat is the file name? I.e Data1.txt : ")

    data_filename = input("\nWhat is the filename you want the data to be stored in called? Please include .txt at the end. I.e 'Data.txt' etc. : ")

    pop_size = int(input("\nWhat is the population size? Please enter an even number: "))
    
    max_gens = int(input("\nWhat is the maximum number of generations the algorithm should run for? Please enter a number: "))

    elitism_size = int(input("\nWhat is the number of 'elite' individuals that should be kept per generation? Please enter a number: "))

    crossover_rate = float(input("\nWhat is the crossover rate? Please enter a float between 0-1 (i.e 0 = 0%, 0.5 = 50%, 1 = 100%): "))

    mutation_rate = float(input("\nWhat is the mutation rate? Please enter a float between 0-1 (i.e 0 = 0%, 0.5 = 50%, 1 = 100%): "))

    crossover_type = input("\nWhat type of crossover should be used? Type in the letter (A/B) \nA. Uniform Order\nB. Two-Point ").lower()

    print_gen = int(input("\nFinally, how often would you like to print out the generations? (i.e every 50, 100, 150 generations, etc.) This will also affect the graph data, as it it will update every X generations. Please enter a number: "))
 

    #data1GA = GeneticAlgorithm("Assignment1/Data1.txt", pop_size, max_gens, elitism_size, crossover_rate, mutation_rate, crossover_type)
    #data2GA = GeneticAlgorithm("Assignment1/Data2.txt", pop_size, max_gens, elitism_size, crossover_rate, mutation_rate, crossover_type)

    GeneticAlgo = GeneticAlgorithm(f"Assignment1\{filename}", pop_size, max_gens, elitism_size, crossover_rate, mutation_rate, crossover_type, print_gen, data_filename)
    
    
    #test = GeneticAlgorithm("Assignment1/Data2.txt", 1, 1, 1, 1, 1.0, 1.0)

