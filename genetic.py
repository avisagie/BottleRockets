
import random
import copy
from operator import itemgetter


verbose = True


def cap(val, lower, upper):
    if lower > upper:
        raise Exception(f'{lower} > {upper}')    
    return max(min(val, upper), lower)


class Param:

    def __init__(self, lower=None, upper=None, param=None):
        """
        Initialize randomly. If param is given, initialize from it. Else use upper and lower.
        """

        if param is not None:

            self.lower = param.lower
            self.upper = param.upper
            self.var = param.var

        else:

            if lower > upper:
                raise Exception(f'{lower} > {upper}')    

            self.lower = lower
            self.upper = upper
            self.var = 0.5 * (upper - lower)

        self.val = cap( (self.upper-self.lower) * random.random() + self.lower, self.lower, self.upper )


    def mutate(self, temperature):
        """
        Perturb the value using the a gaussian with var = temperature * 0.25 * (upper-lower)
        """
        new_val = cap(self.val + random.normalvariate(0.0, temperature * self.var), self.lower, self.upper)
        self.val = new_val


class Minion:

    def __init__(self, params):
        self.params = params


    def randomize(self):
        x = {}
        for k,v in self.params.items():
            if isinstance(v, Param):
                x[k] = Param(param=v)
            else:
                x[k] = v

        self.params = x


    def mutate(self, temperature):
        x = copy.deepcopy(self.params)
        for k, v in x.items():
            if isinstance(v, Param):
                v.mutate(temperature)
        self.params = x


    def crossover(self, other):
        x = {}
        for k in self.params:
            if isinstance(self.params[k], Param):
                x[k] = copy.deepcopy(random.choice([self.params[k], other.params[k]]))
            else:
                x[k] = self.params[k]

        return Minion(x)


    def unpack(self):
        x = {}
        for k in self.params:
            if isinstance(self.params[k], Param):
                x[k] = self.params[k].val
            else:
                x[k] = self.params[k]

        return x


class GeneticAlgorithm:

    def __init__(self, params, fitness, population_size=20, generations=20):
        """
        A genetic algorithm based on parameters params. Params is a dict with 
        string keys and arbitrary values. Only values of type Param will 
        participate in the GA.

        fitness is a function that takes params and return a float. Bigger is 
        better.
        """
        self.params = params
        self.population_size = population_size
        self.generations = generations
        self.population = []
        self.fitness = fitness


    def initialize_population(self):
        if verbose: print("Initializing")

        ret = []
        for ii in range(self.population_size):
            minion = Minion(self.params)
            minion.randomize()
            fitness = self.fitness(minion.unpack())
            ret.append( (fitness, minion) )

        ret.sort(key=itemgetter(0), reverse=True)

        self.population = ret

    def run(self):
        try:
            self.__run()
        except KeyboardInterrupt:
            pass
            
        fitness, minion = self.population[0]
        print()
        print('Stopped.')
        print(f'Best parameters found so far (fitness: {fitness}):')
        for k in sorted(minion.params):
            v = minion.params[k]
            if isinstance(v, Param):
                print(f'    *"{k}": {v.val}')
            else:
                print(f'     "{k}": {v}')


    def __run(self):

        self.initialize_population()

        temperature = 1.0

        for ii in range(self.generations):
            survivors = self.population[:self.population_size // 4]

            if verbose: 
                print(f"Generation {ii+1}, temperature:{temperature}, {len(survivors)} survivors:")
                for fitness, _ in survivors:
                    print(f"    {fitness:0.02f}")

            next_population = []
            for x in survivors:
                next_population.append(x)

            for ii in range(self.population_size // 3):
                _, minion = random.choice(survivors)
                new_minion = copy.deepcopy(minion)
                new_minion.mutate(temperature)
                next_population.append( (self.fitness(new_minion.unpack()), new_minion) )

            while len(next_population) < self.population_size:
                _, mom = random.choice(next_population)
                _, dad = random.choice(next_population)
                baby = mom.crossover(dad)
                next_population.append( (self.fitness(baby.unpack()), baby) )

            next_population.sort(key=itemgetter(0), reverse=True)
            self.population = next_population

            temperature *= 0.9
