#include "BinaryGASolver.h"
#include "simple_svg_1.0.0.hpp"
#include "Common.h"
#include <iostream>
#include <iomanip>
#include <chrono>

// http://www.ai-junkie.com/ga/intro/gat3.html
// Given an area that has a number of non overlapping disks scattered about its surface
// Use a genetic algorithm to find the disk of largest radius which may be placed amongst 
// these disks without overlapping any of them.

// https://code.google.com/archive/p/simple-svg/

#define POPULATION_SIZE 600
#define MUTATION_PROBABILITY 0.01
#define CROSSOVER_FACTOR 0.75

#define MAX_NUMBER_OF_GENERATIONS 15000
#define STOP_AFTER_NUM_GENERATIONS_WITHOUT_CHANGE 1000

#define CIRCLES_SVG_WIDTH 512
#define CIRCLES_SVG_HEIGHT 512
#define CIRCLES_NUMBER 20
#define CIRCLES_RADIUS_FROM 5
#define CIRCLES_RADIUS_TO 50

#define CIRCLES_MUTATION_MAX_VALUE 5.0

struct Circle
{
	Common::Point p;
	double r;
};

struct Problem
{
	std::vector<Circle> circles;
	double width;
	double height;
};

bool ValidateCircle(const Problem & problem, const Circle & circle)
{
	// validate if circle is within frame
	if (circle.p.x < circle.r || circle.p.y < circle.r ||
		circle.p.x > problem.width - circle.r || circle.p.y > problem.height - circle.r)
		return false;
	// validate if circle does not intersect other circle
	for (const auto & problemCircle : problem.circles)
	{
		double distance = Common::Distance(problemCircle.p, circle.p);
		if (distance < problemCircle.r + circle.r)
			return false;
	}
	return true;
}

Problem GenerateProblem(double width, double height, size_t circlesNumber, double radiusFrom, double radiusTo)
{
	Problem problem;
	problem.width = width;
	problem.height = height;

	for (size_t i = 0; i < CIRCLES_NUMBER; ++i)
	{
		while (true)
		{
			Circle circle;
			circle.r = Common::Frand(radiusFrom, radiusTo);
			circle.p.x = Common::Frand(circle.r, width - circle.r);
			circle.p.y = Common::Frand(circle.r, height - circle.r);

			if (ValidateCircle(problem, circle))
			{
				problem.circles.push_back(std::move(circle));
				break;
			}
		}
	}

	return problem;
}

void VisualizeProblem(const char * name, const Problem & problem, Circle * solution = nullptr)
{
	svg::Document doc(name,
		svg::Layout({ problem.width, problem.height }, svg::Layout::BottomLeft));

	for (const auto & circle : problem.circles)
	{
		doc << svg::Circle(svg::Point(circle.p.x, circle.p.y), 2.0*circle.r, 
			svg::Fill(svg::Color(100, 200, 120)));
	}

	if (solution)
	{
		doc << svg::Circle(svg::Point(solution->p.x, solution->p.y), 2.0*solution->r,
			svg::Fill(svg::Color(200, 100, 120)));
	}

	doc << svg::Rectangle({ 0.0, problem.height }, problem.width, problem.height, 
		svg::Fill(), svg::Stroke(1.0, svg::Color::Silver));

	doc.save();
}

Circle ConvertToCircle(const std::vector<double> & chromosome)
{
	Circle ret;

	ret.p.x = chromosome[0];
	ret.p.y = chromosome[1];
	ret.r = chromosome[2];

	return ret;
}

struct EvaluateCircles
{
	BinaryGA::EvaluationResult operator()(uint32_t generation, const std::vector<double> & chromosome)
	{
		if (generation != currentGeneration)
		{
			std::cout << "\rGeneration " << std::fixed << generation;
			std::cout << " max radius " << std::fixed << std::setprecision(2) << maxSolution.r;
			//closestNumber = 0.0;
			currentGeneration = generation;
			numberOfGenerationsWithCurrentSolution++;
		}
		auto circle = ConvertToCircle(chromosome);
		if (ValidateCircle(problem, circle))
		{
			if (circle.r > maxSolution.r)
			{
				maxSolution = circle;
				numberOfGenerationsWithCurrentSolution = 0;

				return BinaryGA::EvaluationResult::ContinueProcessing;
			}
		}

		return numberOfGenerationsWithCurrentSolution >= STOP_AFTER_NUM_GENERATIONS_WITHOUT_CHANGE ?
			BinaryGA::EvaluationResult::ObjectiveReached : BinaryGA::EvaluationResult::ContinueProcessing;
	}
	const Problem & problem;
	uint32_t currentGeneration = 0;
	uint32_t numberOfGenerationsWithCurrentSolution = 0;
	Circle maxSolution;
};

void TestCircles()
{
	std::cout << "Circles" << std::endl;

	auto problem = GenerateProblem(CIRCLES_SVG_WIDTH, CIRCLES_SVG_HEIGHT, 
		CIRCLES_NUMBER, CIRCLES_RADIUS_FROM, CIRCLES_RADIUS_TO);
	VisualizeProblem("data\\CirclesProblem.svg", problem);

	BinaryGA::Definition<double> definition;

	definition.parentSelection = BinaryGA::ParentSelectionType::Ranked;
	definition.mutation = BinaryGA::MutationType::Custom;
	definition.crossover = BinaryGA::CrossoverType::OnePoint;

	definition.populationSize = POPULATION_SIZE;
	definition.mutationProbability = MUTATION_PROBABILITY;
	definition.crossoverFactor = CROSSOVER_FACTOR;
	definition.maxNumberOfGenerations = MAX_NUMBER_OF_GENERATIONS;
	definition.numberOfGenes = 3; // this is number of bits of 3 doubles

	definition.initializationCustomCallback = [](size_t ) -> std::vector<double>
	{
		double radius = Common::Frand(0.0, CIRCLES_SVG_WIDTH);
		return { Common::Frand(radius, CIRCLES_SVG_WIDTH - radius), Common::Frand(radius, CIRCLES_SVG_HEIGHT - radius), radius };
	};

	definition.computeFitness = [&problem](const std::vector<double> & chromosome) -> double
	{
		auto circle = ConvertToCircle(chromosome);
		if (ValidateCircle(problem, circle))
			return circle.r;
		return 0.0; // TODO
	};

	definition.mutationCustomCallback = [](const double & value, size_t index) -> double
	{
		//double upperBound = index == 1 ? CIRCLES_SVG_HEIGHT : CIRCLES_SVG_WIDTH;
		//return value + Common::Frand(0.0, upperBound);
		return value + Common::Frand(-CIRCLES_MUTATION_MAX_VALUE, CIRCLES_MUTATION_MAX_VALUE);
	};

	EvaluateCircles evaluate{ problem };
	definition.evaluate = std::ref(evaluate);

	auto startTime = std::chrono::high_resolution_clock::now();
	auto solution = BinaryGA::Solve(definition);
	std::chrono::duration<double, std::milli> solveDuration = std::chrono::high_resolution_clock::now() - startTime;

	std::cout << std::endl << "Generation " << evaluate.currentGeneration << " (" << solveDuration.count() << "ms)" << std::endl;
	std::cout << "Best found solution: " << std::endl;
	std::cout << "radius: " << std::fixed << std::setprecision(2) << evaluate.maxSolution.r << std::endl;
	std::cout << std::endl;

	VisualizeProblem("data\\CirclesSolution.svg", problem, &evaluate.maxSolution);
}
