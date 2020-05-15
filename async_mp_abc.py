import os
import random
import asyncio
import pandas as pd
import numpy as np

from scipy import stats
from concurrent.futures import ThreadPoolExecutor
from multiple_pulses import WrightFisherPopulation

class AsyncMpABC:

    def __init__(self, workers=22):
        self._executor = ThreadPoolExecutor(workers)


    async def in_thread(self, func, x, populationPulses, m, parentalPopulations, founderPopulation, generations, firstChromosome, lastChromosome,sizePopulationObserved,sizePopulationSimulated, simulatorPath, observedDistribution, output):
        loop = asyncio.get_event_loop()
        return await loop.run_in_executor(self._executor, func, x, populationPulses, m, parentalPopulations, founderPopulation, generations, firstChromosome, lastChromosome,sizePopulationObserved,sizePopulationSimulated, simulatorPath, observedDistribution, output)

    def calculatePopulationAncestry(self, observed_path, parental):
    
        population_ancestry = {}
        
        observed = pd.read_csv(observed_path, sep="\t", header=None)
        observed.columns = ["Ancestry","Size","Individual","Chr"]    
            
        grouped_observed = observed[['Ancestry','Size']].groupby("Ancestry").sum()
        total_size = grouped_observed['Size'].sum() 
        
        for pop in parental:
            population_ancestry[pop] = grouped_observed.loc[pop, 'Size']/total_size

        return population_ancestry
    
    def getGenerations(self, first_generation, no_pulses, parental):
    
        generations_distance = int(first_generation//no_pulses)
        parental_size = len(parental)
        
        total_range = (no_pulses * parental_size + generations_distance + (parental_size-1))
        
        generation_list = [ (x+1) for x in range(total_range-1, 0, -1) if x % generations_distance != 0]
        #generation_list.insert(0, firstGeneration+1)
        
        return generation_list
    
    def getProportion(self, no_pulses, parental):
        
        alpha = [1 for _ in range(len(parental))]
        
        distribution = [np.random.dirichlet(alpha) for _ in range(no_pulses)]
        distribution = np.reshape(distribution, (1,-1))
        distribution = np.append(distribution[0], 1.0)
    
        return distribution
    
    def getPopulation(self, no_pulses, parental, founder):
    
        total_iterations = no_pulses * len(parental) -1
        
        population_pulses = [parental[random.randint(0,len(parental)-1)] for x in range(total_iterations)]
        population_pulses.insert(0, founder)
        
        return population_pulses
    
    def compareDistribution(self, observed, simulated):
    
        D,pvalue=stats.ks_2samp(observed.Size, simulated.Size)
        return pvalue, D

    def executeCommand(self, populationPulses, m, generations, firstChromosome, lastChromosome,sizePopulationObserved, sizePopulationSimulated, simulatorPath, output):

        chromosomes_size = {
            "1": 2.84249773, "2": 2.68831692, "3": 2.23257669,
            "4": 2.14204523, "5": 2.04047482, "6": 1.91826989,
            "7":1.87149372, "8": 1.6800111, "9": 1.66298026,
            "10": 1.80954204, "11": 1.58216679, "12": 1.74513019,
            "13": 1.25475705, "14": 1.18063515, "15": 1.41332565,
            "16": 1.34024191, "17": 1.284838, "18": 1.17473439,
            "19": 1.0773193, "20": 1.08212234, "21": 0.61911307,
            "22": 0.72488198,
        }
        
        for chromosome in range(firstChromosome, lastChromosome+1):
            chromo = chromosome
            r = float(chromosomes_size[str(chromosome)])

            wf = WrightFisherPopulation(population_size=sizePopulationSimulated, migration_times=generations,
                                    migration_probabilities=m,
                                    migration_sources=populationPulses)
            
            sources, end_points = wf.simulate_tracts(r=r)
            tract_lengths = np.diff(end_points)

            response_list = [({"source": i, "tract_lengths": tract_lengths[i], "chromossome": chromo}) for i in range(len(sources))]
        
        return response_list
    
    def execute_sample(self, pulse, populationPulses, m, parentalPopulations, founderPopulation, generations, firstChromosome, lastChromosome,sizePopulationObserved,sizePopulationSimulated, simulatorPath, observedDistribution, output):
        
        print("Executing pulse: {}".format(pulse))

        return self.executeCommand(populationPulses, m, generations, firstChromosome, lastChromosome,sizePopulationObserved,sizePopulationSimulated, simulatorPath, output)
        

    async def concatenate_responses(self, no_pulses, populationPulses, m, parentalPopulations, founderPopulation, generations, firstChromosome, lastChromosome,sizePopulationObserved,sizePopulationSimulated, simulatorPath, observedDistribution, output):
        results = await asyncio.gather( *(self.in_thread(self.execute_sample, x, populationPulses, m, parentalPopulations, founderPopulation, generations, firstChromosome, lastChromosome,sizePopulationObserved,sizePopulationSimulated, simulatorPath, observedDistribution, output) for x in range(50000) ))

        simulatedDistribution=pd.read_csv(results, sep="\t", header=None)
        simulatedDistribution.columns=["Ancestry","Size","Individual","Chr"]
         
        
        #Peform ks test for each ancestry
        response = {}
        #response['pulse'] = pulse
        response['populationPriori'] = populationPulses
        response['mPriori'] = m
        response['pvalue'] = {}
        response['ds'] = {}

        for pop in parentalPopulations:
            observedFiltered = observedDistribution.loc[observedDistribution["Ancestry"] == pop]
            simulatedFiltered = simulatedDistribution.loc[simulatedDistribution["Ancestry"] == pop]

            pvalueTemp,DTemp = self.compareDistribution(observedFiltered, simulatedFiltered)

            if pop in response['pvalue']:
                current_pvalue = response['pvalue']
                current_pvalue = response[pop].append(pvalueTemp)
            else:
                current_pvalue = response['pvalue']
                current_pvalue[pop] = []

            if pop in response['ds']:
                current_pvalue = response['ds']
                current_pvalue = response[pop].append(DTemp)
            else:
                current_pvalue = response['ds']
                current_pvalue[pop] = []
            

        return response

    def run_async_model(self, numberOfPulses, simulatorPath, toSample, firstChromosome, lastChromosome, firstGeneration,sizePopulationObserved, sizePopulationSimulated, founderPopulation, observedFile, parentalPopulations, observedDistribution, output):
        
        loop = asyncio.get_event_loop()
        
        try:

            
            #Observed Population Ancestry
            #M = self.calculatePopulationAncestry(observedFile, parentalPopulations)
    
            #Creating a generation list
            generations = self.getGenerations(firstGeneration, numberOfPulses, parentalPopulations)

            m = self.getProportion(numberOfPulses, parentalPopulations)
            populationPulses = self.getPopulation(numberOfPulses, parentalPopulations, founderPopulation)

            response = loop.run_until_complete(self.concatenate_responses(numberOfPulses, populationPulses, m, parentalPopulations, founderPopulation, generations, firstChromosome, lastChromosome,sizePopulationObserved,sizePopulationSimulated, simulatorPath, observedDistribution, output))
            
            return response
        
        finally:
            loop.run_until_complete(loop.shutdown_asyncgens())
            loop.close()

if __name__ == "__main__":
    parental= ['AFR', 'EUR', 'NAT']
    observed_distribution_path = 'DesenvolvimentoOutSimulatedData'

    observedDistribution = pd.read_csv("Desenvolvimento.txt", sep="\t", header=None)

    AsyncABC = AsyncMpABC(20)
    abc_simulation = AsyncABC.run_async_model(5, "python multiple_pulses.py", 50000, 22, 22, 20, 10000, 100, "NAT", "Desenvolvimento.txt", parental, observedDistribution, "DesenvolvimentoOutSimulatedData")
    
    print(abc_simulation)


'''
--parental AFR NAT EUR --numberOfPulses 5 
#--output DesenvolvimentoOut --input Desenvolvimento.txt --founder NAT --numberOfIndividuals 100 
#--sizePop 10000 --firstGeneration 20 --firstChromosome 22 --lastChromosome 22 --mode 1 --sample 50000 
#--simulator "python multiple_pulses.py" --color blue darkgreen red
'''