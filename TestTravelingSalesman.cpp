#include <fstream>
#include <iostream>
#include <cmath>
#include <map>
#include <chrono>
#include "json.hpp"
#include "BinaryGASolver.h"
#include "Common.h"

using json = nlohmann::json;

#define POPULATION_SIZE 300
#define MUTATION_PROBABILITY 0.01
#define CROSSOVER_FACTOR 0.75

#define MAX_NUMBER_OF_GENERATIONS 15000

using DistancesMap = std::map<std::pair<size_t, size_t>, double>;
double ComputeDistance(const std::vector<uint8_t> & chromosome, const DistancesMap & distances)
{
	double totalDistance = 0.0;
	for (size_t i = 1; i < chromosome.size(); ++i)
		totalDistance += distances.at({ chromosome[i - 1], chromosome[i] });
	totalDistance += distances.at({ chromosome.back(), chromosome[0] });

	return totalDistance;
}

struct EvaluateTravelingSalesman
{
	BinaryGA::EvaluationResult operator()(uint32_t generation, const std::vector<uint8_t> & chromosome)
	{
		if (generation != currentGeneration)
		{
			std::cout << "\rGeneration " << std::fixed << generation;
			std::cout << " current minimum: " << std::fixed << std::setprecision(2) << currentDistance;
			std::cout << " optimal minimum: " << std::fixed << std::setprecision(2) << optimalDistance;
			//currentDistance = std::numeric_limits<double>::max();
			currentGeneration = generation;
		}
		double distance = ComputeDistance(chromosome, distances);
		if (distance < currentDistance)
		{
			currentDistance = distance;
			currentSolution = chromosome;
		}

		return fabs(currentDistance - optimalDistance) < 0.1 ?
			BinaryGA::EvaluationResult::ObjectiveReached : BinaryGA::EvaluationResult::ContinueProcessing;
	}
	std::map<std::pair<size_t, size_t>, double> & distances;
	double optimalDistance;
	uint32_t currentGeneration = 0;
	double currentDistance = std::numeric_limits<double>::max();
	std::vector<uint8_t> currentSolution;
};

std::string ConvertSolutionToString(const std::vector<uint8_t> & solution, const DistancesMap & distances)
{
	std::stringstream str;
	for (size_t i = 0; i < solution.size(); ++i)
		str << (size_t)solution[i] << " -> ";
	str << (size_t)solution[0];
	str << " = " << std::fixed << std::setprecision(3) << ComputeDistance(solution, distances);
	return str.str();
}

void TestTravelingSalesman()
{
	std::cout << "Traveling salesman" << std::endl;

	std::ifstream file("data\\travelingSalesman.json");
	json jsonProblems;
	file >> jsonProblems;

	BinaryGA::Definition<uint8_t> definition;

	definition.parentSelection = BinaryGA::ParentSelectionType::Ranked;
	definition.mutation = BinaryGA::MutationType::Swap;
	definition.crossover = BinaryGA::CrossoverType::Ordered;

	definition.populationSize = POPULATION_SIZE;
	definition.mutationProbability = MUTATION_PROBABILITY;
	definition.crossoverFactor = CROSSOVER_FACTOR;
	definition.maxNumberOfGenerations = MAX_NUMBER_OF_GENERATIONS;

	for (size_t i = 0; i < jsonProblems["problems"].size(); ++i)
	{
		std::cout << "Problem: " << i << std::endl;
		std::cout << "Generation 0";

		// read problem
		double optimalDistance = jsonProblems["problems"][i]["optimal"];
		std::vector<Common::Point> points;
		for (auto point : jsonProblems["problems"][i]["points"])
			points.push_back({ point["x"], point["y"] });

		// precompute distances
		std::map<std::pair<size_t, size_t>, double> distances;
		for (size_t i = 0; i < points.size(); ++i)
		{
			for (size_t j = i + 1; j < points.size(); ++j)
			{
				double distance = Common::Distance(points[i], points[j]);
				distances[{i, j}] = distance;
				distances[{j, i}] = distance;
			}
		}

		// prepare seed
		std::vector<uint8_t> seed;
		for (uint8_t i = 0; i < points.size(); ++i)
			seed.push_back(i);

		definition.initializationCustomCallback = [&seed](size_t) -> std::vector<uint8_t>
		{
			std::random_shuffle(std::begin(seed), std::end(seed));
			return seed;
		};

		definition.numberOfGenes = seed.size();

		// solve

		definition.computeFitness = [&distances](const std::vector<uint8_t> & chromosome) -> double
		{
			return 1.0 / ComputeDistance(chromosome, distances);
		};

		EvaluateTravelingSalesman evaluate{ distances, optimalDistance };
		definition.evaluate = std::ref(evaluate);

		auto startTime = std::chrono::high_resolution_clock::now();
		auto solution = BinaryGA::Solve(definition);
		std::chrono::duration<double, std::milli> solveDuration = std::chrono::high_resolution_clock::now() - startTime;

		std::cout << std::endl << "Generation " << evaluate.currentGeneration << " (" << solveDuration.count() << "ms)" << std::endl;
		if (!solution.empty())
		{
			std::cout << "Optimal solution found: " << std::endl;
			std::cout << ConvertSolutionToString(solution, distances) << std::endl;
		}
		else
		{
			std::cout << "Best found solution: " << std::endl;
			std::cout << ConvertSolutionToString(evaluate.currentSolution, distances) << " ";
			std::cout << std::fixed << std::setprecision(2);
			std::cout << (optimalDistance / (double)ComputeDistance(evaluate.currentSolution, distances)) * 100.0 << "%" << std::endl;
		}
	}

	std::cout << std::endl;
}
