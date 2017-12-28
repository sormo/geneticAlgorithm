#include <fstream>
#include <iostream>
#include <chrono>
#include "json.hpp"
#include "BinaryGASolver.h"

using json = nlohmann::json;

// http://www.ai-junkie.com/ga/intro/gat3.html
// Given the digits 0 through 9 and the operators + , -, * and / , find a sequence 
// that will represent a given target number.The operators will be applied sequentially 
// from left to right as you read.
// Examples:
// target number 23 -> the sequence 6 + 5 * 4 / 2 + 1 would be one possible solution.
// target number 75.5 -> then 5 / 2 + 9 * 7 - 5 would be a possible solution.

#define POPULATION_SIZE 300
#define MUTATION_PROBABILITY 0.01
#define CROSSOVER_FACTOR 0.75

#define MAX_NUMBER_OF_GENERATIONS 15000

struct BinToChar
{
	char bin;
	char letter;
};
static const BinToChar CHARACTER_MAP[] =
{
	{ 0b0000, '0' },
	{ 0b0001, '1' },
	{ 0b0010, '2' },
	{ 0b0011, '3' },
	{ 0b0100, '4' },
	{ 0b0101, '5' },
	{ 0b0110, '6' },
	{ 0b0111, '7' },
	{ 0b1000, '8' },
	{ 0b1001, '9' },
	{ 0b1010, '+' },
	{ 0b1011, '-' },
	{ 0b1100, '*' },
	{ 0b1101, '/' },
};
#define COUNTOF(c) (sizeof(c)/(sizeof(c[0])))
#define CHARACTER_BIN_SIZE 4
#define MAX_NUMBER_OF_CHARACTERS 40 // a + b + c + d ...
#define IS_NUMBER(c) (c < 10)

char ConvertBinCharToChar(char binChar)
{
	for (size_t i = 0; i < COUNTOF(CHARACTER_MAP); ++i)
	{
		if (CHARACTER_MAP[i].bin == binChar)
			return CHARACTER_MAP[i].letter;
	}
	return 'X';;
}

void ApplyOperator(char op, char value, uint32_t & out)
{
	switch (op)
	{
	case 0b1010: // +
		out += value;
		break;
	case 0b1011: // -
		out -= value;
		break;
	case 0b1100: // *
		out *= value;
		break;
	case 0b1101: // /
		if (value != 0)
			out /= value;
		break;
	default:
		out = value;
	}
}

uint32_t ConvertToNumber(const std::vector<bool> & binary, std::function<void(char)> cbk = nullptr)
{
	uint32_t ret = 0;
	bool needNumber = true;
	char lastOperator = 0;

	for (size_t i = 0; i < binary.size(); i += CHARACTER_BIN_SIZE)
	{
		char binChar = 0;
		for (size_t j = 0; j < CHARACTER_BIN_SIZE; ++j)
		{
			if (binary[i+j])
				binChar |= 1 << (CHARACTER_BIN_SIZE - 1 - j);
		}

		if (binChar >= COUNTOF(CHARACTER_MAP))
			continue;

		if ((needNumber && IS_NUMBER(binChar)) || (!needNumber && !IS_NUMBER(binChar)))
		{
			if (needNumber)
				ApplyOperator(lastOperator, binChar, ret);
			else
				lastOperator = binChar;
			needNumber = !needNumber;

			if (cbk)
				cbk(binChar);
		}
	}

	return ret;
}

std::string ConvertToString(const std::vector<bool> & binary)
{
	std::string ret;
	char lastOperator = 0;
	ConvertToNumber(binary, [&ret, &lastOperator](char binChar)
	{
		if (IS_NUMBER(binChar))
		{
			if (lastOperator)
			{
				ret.push_back(ConvertBinCharToChar(lastOperator));
				ret.push_back(' ');
			}
			ret.push_back(ConvertBinCharToChar(binChar));
			ret.push_back(' ');
		}
		else
			lastOperator = binChar;
	});
	return ret;
}

struct EvaluateDigitsAndOperators
{
	BinaryGA::EvaluationResult operator()(uint32_t generation, const std::vector<bool> & chromosome)
	{
		if (generation != currentGeneration)
		{
			std::cout << "\rGeneration " << std::fixed << generation;
			std::cout << " value " << std::fixed << std::setprecision(2) << (closestNumber / targetNumber) * 100.0;
			//closestNumber = 0.0;
			currentGeneration = generation;
		}
		double number = ConvertToNumber(chromosome);
		if (fabs(number - targetNumber) < fabs(closestNumber - targetNumber))
		{
			closestNumber = number;
			closestSolution = chromosome;
		}

		return fabs(number - targetNumber) < 0.0001 ?
			BinaryGA::EvaluationResult::ObjectiveReached : BinaryGA::EvaluationResult::ContinueProcessing;
	}
	uint32_t currentGeneration = 0;
	double closestNumber = 0.0;
	std::vector<bool> closestSolution;
	double targetNumber = 0.0;
};

void TestDigitsAndOperators()
{
	std::cout << "Digits and operators" << std::endl;

	std::ifstream file("data\\digitsAndOperators.json");
	json jsonProblems;
	file >> jsonProblems;

	BinaryGA::Definition<bool> definition;

	definition.parentSelection = BinaryGA::ParentSelectionType::Ranked;
	definition.mutation = BinaryGA::MutationType::Toggle;
	definition.crossover = BinaryGA::CrossoverType::OnePoint;

	definition.populationSize = POPULATION_SIZE;
	definition.mutationProbability = MUTATION_PROBABILITY;
	definition.crossoverFactor = CROSSOVER_FACTOR;
	definition.maxNumberOfGenerations = MAX_NUMBER_OF_GENERATIONS;
	definition.numberOfGenes = MAX_NUMBER_OF_CHARACTERS * CHARACTER_BIN_SIZE;

	for (size_t i = 0; i < jsonProblems["problems"].size(); ++i)
	{
		std::cout << "Problem: " << i << std::endl;
		std::cout << "Generation 0";

		// read problem
		double targetNumber = jsonProblems["problems"][i];

		definition.computeFitness = [targetNumber](const std::vector<bool> & chromosome) -> double
		{
			return 1.0 / (targetNumber - (double)ConvertToNumber(chromosome));
		};

		EvaluateDigitsAndOperators evaluate;
		evaluate.targetNumber = targetNumber;
		definition.evaluate = std::ref(evaluate);

		auto startTime = std::chrono::high_resolution_clock::now();
		auto solution = BinaryGA::Solve(definition);
		std::chrono::duration<double, std::milli> solveDuration = std::chrono::high_resolution_clock::now() - startTime;

		std::cout << std::endl << "Generation " << evaluate.currentGeneration << " (" << solveDuration.count() << "ms)" << std::endl;
		if (!solution.empty())
		{
			std::cout << "Optimal solution found: " << std::endl;
			std::cout << ConvertToString(solution) << "= " << targetNumber << std::endl;
		}
		else
		{
			std::cout << "Best found solution: " << std::endl;
			std::cout << ConvertToString(evaluate.closestSolution) << "= " << targetNumber << " ";
			std::cout << std::fixed << std::setprecision(2) << (evaluate.closestNumber / evaluate.targetNumber) * 100.0 << std::endl;
		}
	}

	std::cout << std::endl;
}
