import sys
import argparse as ap
import numpy as np


def getchar():
    print("Aperte ENTER")
    sys.stdin.read(1)


class WrightFisherPopulation(object):
    def __init__(self, population_size, migration_times,
                 migration_probabilities, migration_sources):
        #       #Tamanho da populacao
        self.population_size = population_size

#	#Verificar se o ultimo local do vetor é igual 1 (founder)
        assert migration_probabilities[-1] == 1.0, (f"Final migration probability must be 1.0. It is {migration_probabilities[-1]}")
#	#Cria dicionarios migration time -> migration probability
        self.migration_probabilities = dict(zip(migration_times,
                                                migration_probabilities))
#	#Cria dicionarios migration time -> migration sources
        self.migration_sources = dict(zip(migration_times, migration_sources))

        self.population = {}

#	#Limite de recursão 3*migration time +100
        sys.setrecursionlimit(3 * migration_times[-1] + 100)

    def simulate_mosaic(self, recombination_distance):
        #	#ancestrais locais = get_person da geracao =1
        local_ancestors = [self.get_person(generation=1)]
#	Enquanto a ultima populacao local é None, isso é, nao é de um migrante
        while local_ancestors[-1].source_population is None:
            #	    adiciona um ancestral através de um get_parent
            local_ancestors.append(local_ancestors[-1].get_parent(0, self))

        # print(local_ancestors)
        tile_sources = []
        tile_boundaries = [0]

        position = local_ancestors[-1].get_tile_length()
#
        while position < recombination_distance:
            # print("position"+str(position))
            #	Escolhe aleatoriamente um dos ancestrais para ser o recombiner
            recombiner = np.random.choice(local_ancestors[:-1])

#	en1quanto o ultimo nao for o recombiner
            while local_ancestors[-1] is not recombiner:
                local_ancestors.pop().last_seen_at = position

#	Se nao eh o recombiner, marco o last_seen_at

# atualiza last seen at e swapa o copying
            recombiner.recombine(position)

# Adiciona individuos ate que o source seja diferente de None
            while local_ancestors[-1].source_population is None:
                ancestor = local_ancestors[-1].get_parent(position, self)
                local_ancestors.append(ancestor)

# Fonte recebe o ancestor
            tile_sources.append(ancestor.source_population)
# Apend na posicao dos boundaries
            tile_boundaries.append(position)
            position += local_ancestors[-1].get_tile_length()
            #print("A source foi"+ancestor.source_population +" e a position eh "+str(position))

# Adiciona o ultimo e boundaries
        tile_sources.append(local_ancestors[-1].source_population)
        tile_boundaries.append(recombination_distance)
        return tile_sources, tile_boundaries

    def simulate_tracts(self, r):
        tile_sources, tile_boundaries = self.simulate_mosaic(r)
        assert len(tile_sources) + 1 == len(tile_boundaries)

        current_source = tile_sources[0]

        tract_sources = [current_source]
        tract_boundaries = [0]
        for source, boundary in zip(tile_sources, tile_boundaries[:-1]):
            if source == current_source:
                continue
            tract_sources.append(source)
            tract_boundaries.append(boundary)
            current_source = source

        tract_boundaries.append(tile_boundaries[-1])

        return tract_sources, tract_boundaries

    def get_person(self, generation):
        #	#person_id= (geração, numero aleatorio até 2*tamanho da populacao)
        person_id = (generation, np.random.randint(2 * self.population_size))

#	Se o id já existe na lista de populacoes
        if person_id in self.population:
            #	    Retorna o individuo
            return self.population[person_id]

#	Caso nao exista
        source_population = None
#	Se a geracao existe no dicionario geracao-> probabilidade
        if generation in self.migration_probabilities:
            #	    Caso seja uma geração com migração, gera um numero aleatorio e se for menor o antepassado do individuo será um que chegou
            if np.random.uniform() < self.migration_probabilities[generation]:
                source_population = self.migration_sources[generation]
#	Se da um return em um indivíduo composto por geracao e source que é none ou migração.
        self.population[person_id] = Person(generation, source_population)
        return self.population[person_id]


class Person(object):
    def __init__(self, generation, source_population):
        #	Geracao
        self.generation = generation
#	Populacao fonte (None ou migrante)
        self.source_population = source_population

#	Insercao de ultima vez visto em -infinito, nao copiou e nao copiando
        self.last_seen_at = -np.inf
        self.copying = None
        self.not_copying = None

    def printObj(self):
        print("Generation: "+str(self.generation))
        print("Source("+str(self.generation)+"): "+str(self.source_population))
        print("Last:("+str(self.generation)+"): "+str(self.last_seen_at))

        if self.copying == None:
            print("Copying ("+str(self.generation)+"): None")
        else:
            print("Copying("+str(self.generation)+"): " +
                  str(self.copying.printObj()))

        if self.not_copying == None:
            print("!Copying:("+str(self.generation)+"): None")
        else:
            print("!Copying: ("+str(self.generation)+"): " +
                  str(self.not_copying.printObj()))

# Swap coping e nao copying
    def swap_copying(self):
        self.copying, self.not_copying = self.not_copying, self.copying

# atualiza last seen at e swapa o copying
    def recombine(self, position):
        self.last_seen_at = position
        self.swap_copying()

    def get_parent(self, position, population):
        distance_from_last_seen = position - self.last_seen_at
        self.last_seen_at = position

#	a probabilidade de swap é 0.5-0.5 * e^(-2*distancia)
#	se o last_seen_at = -np.inf a probabilidade será de -inf
        swap_probability = .5 - .5 * np.exp(-2 * distance_from_last_seen)
        if np.random.uniform() < swap_probability:
            #           Faz o swap copying
            self.swap_copying()
        if self.copying is None:
            #	    Se o copying eh none, get_person
            self.copying = population.get_person(self.generation + 1)
        return self.copying

    def get_tile_length(self):
        assert self.source_population is not None
# Se a source population nao eh none, ele amostra o tamnho de uma exponencial cujo 1/beta = 1/geracaoDeChegada
# geracoes antigas, fragmentos menores. geracoes recentes, fragmentos maiores
        return np.random.exponential(1.0 / self.generation)


if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Simulate migrant tracts.')
    parser.add_argument(
        '-N', help='effective number of diploid people in population')
    parser.add_argument('-r', help='recombination distance, in Morgans')
    parser.add_argument('-m', nargs='+', help='migration probabilities')
    parser.add_argument('-T', nargs='+', help='migration times')
    parser.add_argument('-s', nargs='+', help='source population labels')
    parser.add_argument('-q', help='quantidade de cromossomos', default=1)
    parser.add_argument('-c', help='numero do chr', default=1)
    args = parser.parse_args()
    quantidade = int(args.q)
    cromossomo = str(args.c)
    assert cromossomo != 'X'

    for j in range(quantidade):
        N = int(args.N)
        r = float(args.r)
        sources = args.s
        Ts = np.array(args.T, dtype='int')
        ms = np.array(args.m, dtype='float')

        wf = WrightFisherPopulation(population_size=N, migration_times=Ts,
                                    migration_probabilities=ms,
                                    migration_sources=sources)
        sim_count = 1
        sources, end_points = wf.simulate_tracts(r=r)
        tract_lengths = np.diff(end_points)

        for i in range(len(sources)):
            print ('%s\t%f\t%s\t%s' %
                   (sources[i], tract_lengths[i], j, cromossomo))

