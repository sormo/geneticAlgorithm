#include <fstream>
#include <iostream>
#include <chrono>
#include "BinaryGASolver.h"
#include "json.hpp"

using json = nlohmann::json;

#define POPULATION_SIZE 300
#define MUTATION_PROBABILITY 0.1
#define CROSSOVER_FACTOR 0.75

#define MAX_NUMBER_OF_GENERATIONS 15000

struct Problem
{
	std::vector<uint32_t> profits;
	std::vector<uint32_t> weights;
	std::vector<bool> optimal;
	uint32_t capacity;
};

uint32_t ComputeValue(const std::vector<bool> & flags, const std::vector<uint32_t> & values)
{
	uint32_t ret = 0;
	for (size_t i = 0; i < flags.size(); ++i)
	{
		if (flags[i])
			ret += values[i];
	}
	return ret;
}

std::tuple<uint32_t, uint32_t> ComputeProfitAndWeight(const std::vector<bool> & flags, const Problem & problem)
{
	uint32_t profit = 0;
	uint32_t weight = 0;

	for (size_t i = 0; i < flags.size(); ++i)
	{
		if (flags[i])
		{
			if (weight + problem.weights[i] >= problem.capacity)
				break;
			profit += problem.profits[i];
			weight += problem.weights[i];
		}
	}

	return { profit, weight };
}

struct EvaluateKnapsack
{
	BinaryGA::EvaluationResult operator()(uint32_t generation, const std::vector<bool> & chromosome)
	{
		if (generation != currentGeneration)
		{
			std::cout << "\rGeneration " << std::fixed << generation;
			std::cout << " current profit: " << std::fixed << std::setprecision(2) << currentProfit;
			std::cout << " optimal profit: " << std::fixed << std::setprecision(2) << totalProfit;
			//currentProfit = 0.0;
			currentGeneration = generation;
		}
		double profit = (double)std::get<0>(ComputeProfitAndWeight(chromosome, problem));
		if (profit > currentProfit)
		{
			currentProfit = profit;
			currentSolution = chromosome;
		}

		return problem.optimal == chromosome ?
			BinaryGA::EvaluationResult::ObjectiveReached : BinaryGA::EvaluationResult::ContinueProcessing;
	}
	const Problem & problem;
	double totalProfit;
	double currentProfit = 0.0;
	std::vector<bool> currentSolution;
	uint32_t currentGeneration = 0;
};

std::string ConvertSolutionToString(const std::vector<bool> & solution, const Problem & problem)
{
	std::stringstream str;
	bool first = true;
	for (size_t i = 0; i < solution.size(); ++i)
	{
		if (solution[i])
		{
			if (first)
				first = false;
			else
				str << "+ ";
			str << problem.profits[i] << " (" << problem.weights[i] << ") ";
		}
	}
	auto profitAndWeight = ComputeProfitAndWeight(solution, problem);
	str << "= " << std::get<0>(profitAndWeight) << " (" << std::get<1>(profitAndWeight) << ")";
	return str.str();
}

void TestKnapsack()
{
	std::cout << "Knapsack" << std::endl;

	std::ifstream file("data\\knapsack.json");
	json jsonProblems;
	file >> jsonProblems;

	BinaryGA::Definition<bool> definition;

	definition.populationSize = POPULATION_SIZE;
	definition.mutationProbability = MUTATION_PROBABILITY;
	definition.crossoverFactor = CROSSOVER_FACTOR;
	definition.maxNumberOfGenerations = MAX_NUMBER_OF_GENERATIONS;

	definition.parentSelection = BinaryGA::ParentSelectionType::Ranked;
	definition.mutation = BinaryGA::MutationType::Toggle;
	definition.crossover = BinaryGA::CrossoverType::OnePoint;

	for (size_t i = 0; i < jsonProblems["problems"].size(); ++i)
	{
		std::cout << "Problem: " << i << std::endl;
		std::cout << "Generation 0";

		// read problem
		Problem problem =
		{
			jsonProblems["problems"][i]["profits"],
			jsonProblems["problems"][i]["weights"],
			jsonProblems["problems"][i]["optimal"],
			jsonProblems["problems"][i]["capacity"]
		};

		assert(problem.weights.size() == problem.profits.size());
		
		definition.numberOfGenes = problem.weights.size();
		definition.computeFitness = [&problem](const std::vector<bool> & chromosome) -> double
		{
			return (double)std::get<0>(ComputeProfitAndWeight(chromosome, problem));
		};

		EvaluateKnapsack evaluate{ problem, (double)std::get<0>(ComputeProfitAndWeight(problem.optimal, problem)) };
		definition.evaluate = std::ref(evaluate);

		auto startTime = std::chrono::high_resolution_clock::now();
		auto solution = BinaryGA::Solve(definition);
		std::chrono::duration<double, std::milli> solveDuration = std::chrono::high_resolution_clock::now() - startTime;

		std::cout << std::endl << "Generation " << evaluate.currentGeneration << " (" << solveDuration.count() << "ms)" << std::endl;
		if (!solution.empty())
		{
			std::cout << "Optimal solution found: " << std::endl;
			std::cout << ConvertSolutionToString(solution, problem) << std::endl;
		}
		else
		{
			std::cout << "Best found solution: " << std::endl;
			std::cout << ConvertSolutionToString(evaluate.currentSolution, problem) << std::endl;
			std::cout << std::fixed << std::setprecision(2);
			std::cout << ((double)std::get<0>(ComputeProfitAndWeight(evaluate.currentSolution, problem)) / 
				(double)std::get<0>(ComputeProfitAndWeight(problem.optimal, problem))) * 100.0 << "%" << std::endl;
		}
	}

	std::cout << std::endl;
}
