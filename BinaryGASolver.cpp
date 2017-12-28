#include "BinaryGASolver.h"
#include <algorithm>
#include <unordered_set>
#include "Common.h"

namespace BinaryGA
{
	template<typename T>
	struct Chromosome
	{
		std::vector<T> genes;
		uint32_t age = 0;
		double fitness = 0.0;
	};

	template<typename T>
	using Population = std::vector<Chromosome<T>>;

	template<typename T>
	void Mutate(Chromosome<T> & chromosome, double probability, MutationType type, const CustomMutation<T> & custom);

	template<typename T>
	Population<T> InitializePopulation(const Definition<T> & definition)
	{
		std::vector<Chromosome<T>> population;

		for (size_t i = 0; i < definition.populationSize; ++i)
		{
			Chromosome<T> chromosome;
			if (definition.initializationCustomCallback)
			{
				chromosome.genes = definition.initializationCustomCallback(i);
			}
			else
			{
				chromosome.genes.resize(definition.numberOfGenes);
				Mutate(chromosome, 0.5, definition.mutation, definition.mutationCustomCallback);
			}
			population.push_back(std::move(chromosome));
		}

		return population;
	}

	template<typename T>
	bool CheckTerminationCondition(uint32_t generationNumber, const Definition<T> & definition, const Population<T> & population, std::vector<T> & ret)
	{
		for (const auto & chromosome : population)
		{
			if (definition.evaluate(generationNumber, chromosome.genes) == EvaluationResult::ObjectiveReached)
			{
				ret = chromosome.genes;
				return true;
			}
		}

		if (generationNumber >= definition.maxNumberOfGenerations)
			return true;

		return false;
	}

	// SELECTION //////////////////////////////////////////////////////////////

	template<class T>
	std::vector<const T*> Select(const std::vector<T> & population, size_t number, std::function<double(const T&)> fitness)
	{
		// http://www.obitko.com/tutorials/genetic-algorithms/selection.php
		// [Sum] Calculate sum of all chromosome fitnesses in population - sum S.
		// [Select] Generate random number from interval(0, S) - r.
		// [Loop] Go through the population and sum fitnesses from 0 - sum s.
		//        When the sum s is greater then r, stop and return the chromosome where you are.

		double sumOfFitnesses = 0.0;
		for (const T & chromosome : population)
			sumOfFitnesses += fitness(chromosome);

		std::vector<const T*> ret;

		for (size_t i = 0; i < number; ++i)
		{
			double partialSum = Common::Frand(0.0, sumOfFitnesses);
			double fitnessSum = 0.0;
			for (const T & chromosome : population)
			{
				fitnessSum += fitness(chromosome);
				if (fitnessSum >= partialSum)
				{
					ret.push_back(&chromosome);
					break;
				}
			}
		}

		return ret;
	}

	template<typename T>
	std::vector<const Chromosome<T>*> RouletteWheelSelection(const std::vector<Chromosome<T>> & population, size_t number)
	{
		return Select<Chromosome<T>>(population, number, [](const Chromosome<T> & c) { return c.fitness; });
	}

	template<typename T>
	std::vector<const Chromosome<T>*> RankSelection(const std::vector<Chromosome<T>> & population, size_t number)
	{
		// rank
		std::vector<std::pair<size_t, const Chromosome<T>*>> rankedPopulation;
		for (const auto & chromosome : population)
			rankedPopulation.push_back({ 0, &chromosome });
		std::sort(std::begin(rankedPopulation), std::end(rankedPopulation),
			[](const auto & a, const auto & b) { return a < b; });
		for (size_t i = 0; i < rankedPopulation.size(); ++i)
			rankedPopulation[i].first = i + 1;
		// shuffle
		std::random_shuffle(std::begin(rankedPopulation), std::end(rankedPopulation));

		// TODO sum of fitnesses is equal to sum of first population.size() numbers
		std::vector<const std::pair<size_t, const Chromosome<T>*>*> selectedChromosomes;
		selectedChromosomes = Select<std::pair<size_t, const Chromosome<T>*>>(rankedPopulation, number,
			[](const std::pair<size_t, const Chromosome<T>*> & c) { return c.first; });

		std::vector<const Chromosome<T>*> ret;
		for (const auto & chromosome : selectedChromosomes)
			ret.push_back(chromosome->second);

		return ret;

	}

	template<typename T>
	std::vector<const Chromosome<T>*> SelectParent(size_t crossoverSize, Population<T> & population, const Definition<T> & definition)
	{
		switch (definition.parentSelection)
		{
		case ParentSelectionType::Ranked:
			return RankSelection(population, crossoverSize);
		case ParentSelectionType::RouletteWheel:
			return RouletteWheelSelection(population, crossoverSize);
		}

		return std::vector<const Chromosome<T>*>();
	}

	///////////////////////////////////////////////////////////////////////////

	// CROSSOVER //////////////////////////////////////////////////////////////

	template<typename T>
	std::vector<Chromosome<T>> OnePointCrossover(const Chromosome<T> & first, const Chromosome<T> & second)
	{
		const size_t numberOfGenes = first.genes.size();
		size_t point = rand() % numberOfGenes;

		std::vector<Chromosome<T>> ret(2);
		for (size_t i = 0; i < numberOfGenes; ++i)
		{
			ret[0].genes.push_back(i < point ? first.genes[i] : second.genes[i]);
			ret[1].genes.push_back(i < point ? second.genes[i] : first.genes[i]);
		}

		return ret;
	}

	template<typename T>
	Chromosome<T> OrderedCrossoverCreateChild(size_t point, const Chromosome<T> & first, const Chromosome<T> & second)
	{
		const size_t numberOfGenes = first.genes.size();
		std::unordered_set<T> inserted;
		
		Chromosome<T> ret;
		for (size_t i = 0; i < point; ++i)
		{
			ret.genes.push_back(first.genes[i]);
			inserted.insert(first.genes[i]);
		}

		for (size_t i = 0; i < numberOfGenes; ++i)
		{
			if (inserted.find(second.genes[i]) == std::end(inserted))
				ret.genes.push_back(second.genes[i]);
		}

		return ret;
	}

	template<typename T>
	std::vector<Chromosome<T>> OrderedCrossover(const Chromosome<T> & first, const Chromosome<T> & second)
	{
		const size_t numberOfGenes = first.genes.size();
		size_t point = rand() % numberOfGenes;
		
		std::vector<Chromosome<T>> ret{ 
			OrderedCrossoverCreateChild(point, first, second), 
			OrderedCrossoverCreateChild(point, first, second) 
		};

		return ret;
	}

	template<typename T>
	std::vector<Chromosome<T>> Crossover(const Chromosome<T> & first, const Chromosome<T> & second, const Definition<T> & definition)
	{
		if (definition.crossover == CrossoverType::None)
			return { first, second };

		switch (definition.crossover)
		{
		case CrossoverType::OnePoint:
			return OnePointCrossover(first, second);
		case CrossoverType::Ordered:
			return OrderedCrossover(first, second);
		}

		return std::vector<Chromosome<T>>();
	}

	///////////////////////////////////////////////////////////////////////////

	// MUTATION ///////////////////////////////////////////////////////////////

	template<typename T>
	void MutationToggle(Chromosome<T> & chromosome, double mutationProbability)
	{
		for (size_t i = 0; i < chromosome.genes.size(); ++i)
		{
			double mutate = (double)rand() / (double)RAND_MAX;
			if (mutate < mutationProbability)
				chromosome.genes[i] = !chromosome.genes[i];
		}
	}

	template<typename T>
	void MutationSwap(Chromosome<T> & chromosome, double mutationProbability)
	{
		for (size_t i = 0; i < chromosome.genes.size(); ++i)
		{
			double mutate = (double)rand() / (double)RAND_MAX;
			if (mutate < mutationProbability)
			{
				size_t swapIndex = rand() % chromosome.genes.size();
				std::swap(*(std::begin(chromosome.genes) + i), *(std::begin(chromosome.genes) + swapIndex));
			}
		}
	}

	template<typename T>
	void MutationCustom(Chromosome<T> & chromosome, double mutationProbability, const CustomMutation<T> & custom)
	{
		for (size_t i = 0; i < chromosome.genes.size(); ++i)
		{
			double mutate = (double)rand() / (double)RAND_MAX;
			if (mutate < mutationProbability)
			{
				chromosome.genes[i] = custom(chromosome.genes[i], i);
			}
		}
	}

	template<typename T>
	void Mutate(Chromosome<T> & chromosome, double probability, MutationType type, const CustomMutation<T> & custom)
	{
		switch (type)
		{
		case MutationType::Toggle:
			MutationToggle(chromosome, probability);
			break;
		case MutationType::Swap:
			MutationSwap(chromosome, probability);
			break;
		case MutationType::Custom:
			MutationCustom(chromosome, probability, custom);
			break;
		}
	}

	template<typename T>
	void Mutation(std::vector<Chromosome<T>> & childs, const Definition<T> & definition)
	{
		if (definition.mutation == MutationType::None)
			return;

		for (Chromosome<T> & c : childs)
			Mutate(c, definition.mutationProbability, definition.mutation, definition.mutationCustomCallback);
	}

	///////////////////////////////////////////////////////////////////////////

	template<typename T>
	void SurvivorSelection(Population<T> & population, uint32_t populationSize, std::vector<Chromosome<T>> && childs)
	{
		// add childs to population
		std::move(std::begin(childs), std::end(childs), std::back_inserter(population));

		if (population.size() > populationSize)
		{
			// if population is too big, sort it with ration fitness over age
			std::sort(std::begin(population), std::end(population),
				[](const Chromosome<T> & a, const Chromosome<T> & b)
			{
				double ac = a.age == 0 ? (double)a.fitness : (double)a.fitness / (double)a.age;
				double bc = b.age == 0 ? (double)b.fitness : (double)b.fitness / (double)b.age;
				return ac < bc;
			});

			// remove unlucky ones
			size_t toRemove = population.size() - populationSize;
			for (size_t i = 0; i < toRemove; ++i)
				population.erase(population.begin());
		}
	}

	template<typename T>
	void ValidateDefinition(const Definition<T> & definition)
	{
		if (definition.mutation == MutationType::Custom)
		{
			if (!definition.mutationCustomCallback)
				throw std::runtime_error("custom mutation callback not provided");
		}
	}

	template<typename T>
	std::vector<T> Solve(const Definition<T> & definition)
	{
		ValidateDefinition(definition);

		Population<T> population = InitializePopulation(definition);

		std::vector<T> ret;
		uint32_t generationNumber = 0;
		while (!CheckTerminationCondition(generationNumber, definition, population, ret))
		{
			// increment age
			std::for_each(std::begin(population), std::end(population), [](Chromosome<T> & c) { c.age++; });

			// recompute 
			std::for_each(std::begin(population), std::end(population), 
				[&definition](Chromosome<T> & c) { c.fitness = definition.computeFitness(c.genes); });

			// number of crossovers is computed with crossover factor
			size_t crossoverSize = (size_t)((double)population.size() * definition.crossoverFactor);
			crossoverSize = crossoverSize % 2 ? crossoverSize - 1 : crossoverSize;

			// select parents
			std::vector<const Chromosome<T>*> parents = SelectParent(crossoverSize, population, definition);

			// iterate over pairs of parents
			std::vector<Chromosome<T>> childs;
			for (size_t i = 0; i < parents.size(); i = i + 2)
			{
				// crossover - generate childs
				std::vector<Chromosome<T>> tmpChilds = Crossover(*parents[i], *parents[i + 1], definition);
				// mutation
				Mutation(tmpChilds, definition);
				// evaluate childs
				std::for_each(std::begin(tmpChilds), std::end(tmpChilds),
					[&definition](Chromosome<T> & c) { c.fitness = definition.computeFitness(c.genes); });

				for (auto && c : tmpChilds)
					childs.push_back(std::move(c));
			}
			// survivor selection - merge parents childs - keep constant population size
			SurvivorSelection(population, definition.populationSize, std::move(childs));

			generationNumber++;
		}

		return ret;
	}
}
