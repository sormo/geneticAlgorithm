#include "BinaryGASolver.h"
#include "simple_svg_1.0.0.hpp"
#include "Common.h"
#include <iostream>
#include <iomanip>
#include <chrono>

// find path from bottom of screen up with least nodes
// path must not collide with rectangle

#define POPULATION_SIZE 400
#define MUTATION_PROBABILITY 0.01
#define CROSSOVER_FACTOR 0.75

#define MAX_NUMBER_OF_GENERATIONS 15000
#define STOP_AFTER_NUM_GENERATIONS_WITHOUT_CHANGE 300

#define RECTANGLES_SVG_WIDTH 512
#define RECTANGLES_SVG_HEIGHT 512
#define RECTANGLES_NUMBER 100
#define RECTANGLES_SIZE_FROM 15
#define RECTANGLES_SIZE_TO 60
#define BORDER_ON_BOTTOM_AND_UPPER_PART 50

#define START_POINT { RECTANGLES_SVG_WIDTH/2.0, 0.0 }
#define STEP_MAX_LENGTH 10.0

// maximum number of nodes on path
#define MAXIMUM_PATH_SIZE 150

struct Rectangle
{
	Common::Point upperLeft;
	double width;
	double height;
};

struct Problem
{
	std::vector<Rectangle> rectangles;
	double width;
	double height;
};

using Path = std::vector<Common::Point>;

// https://stackoverflow.com/questions/99353/how-to-test-if-a-line-segment-intersects-an-axis-aligned-rectange-in-2d
bool ValidateLineSegmentAgaintRectangle(const Rectangle & rect, const Common::Point & p1, const Common::Point & p2)
{
	auto ImplicitLineEquationThroughTwoPoints = [&p1, &p2](double x, double y)
	{
		// F(x y) = (y2-y1)*x + (x1-x2)*y + (x2*y1-x1*y2)
		return (p2.y - p1.y)*x + (p1.x - p2.x)*y + (p2.x*p1.y - p1.x * p2.y);
	};

	double upperLeftSide = ImplicitLineEquationThroughTwoPoints(rect.upperLeft.x, rect.upperLeft.y);
	double upperRightSide = ImplicitLineEquationThroughTwoPoints(rect.upperLeft.x + rect.width, rect.upperLeft.y);
	double bottomRightSide = ImplicitLineEquationThroughTwoPoints(rect.upperLeft.x + rect.width, rect.upperLeft.y - rect.height);
	double bottomLeftSide = ImplicitLineEquationThroughTwoPoints(rect.upperLeft.x, rect.upperLeft.y - rect.height);

	// Check if all four corners of the rectangle are on the same side of the line. The implicit equation for a line through p1 and p2 is:
	if ((upperLeftSide > 0.0 && upperRightSide > 0.0 && bottomRightSide > 0.0 && bottomLeftSide > 0.0) ||
		(upperLeftSide < 0.0 && upperRightSide < 0.0 && bottomRightSide < 0.0 && bottomLeftSide < 0.0))
		return true;
	
	Common::Point bottomLeft{ rect.upperLeft.x, rect.upperLeft.y - rect.height };
	Common::Point upperRight{ rect.upperLeft.x + rect.width, rect.upperLeft.y };

	// Project the endpoint onto the x axis, and check if the segment's shadow intersects the polygon's shadow. Repeat on the y axis:
	if (p1.x > upperRight.x && p2.x > upperRight.x)
		return true;
	if (p1.x < bottomLeft.x && p2.x < bottomLeft.x)
		return true;
	if (p1.y > upperRight.y && p2.y > upperRight.y)
		return true;
	if (p1.y < bottomLeft.y && p2.y < bottomLeft.y)
		return true;

	return false;
}

// return true of segment does not intercept any rectangle
// otherwise return false
bool ValidateLineSegment(const Problem & problem, const Common::Point & p1, const Common::Point & p2)
{
	if (p1.x < 0.0 || p1.y < 0.0 || p2.x < 0.0 || p2.y < 0.0)
		return false;
	if (p1.x > problem.width || p1.y > problem.height || p2.x > problem.width || p2.y > problem.height)
		return false;

	for (const auto & rect : problem.rectangles)
	{
		if (!ValidateLineSegmentAgaintRectangle(rect, p1, p2))
			return false;
	}

	return true;
}

Problem GenerateProblem()
{
	Problem problem;
	problem.width = RECTANGLES_SVG_WIDTH;
	problem.height = RECTANGLES_SVG_HEIGHT;

	for (size_t i = 0; i < RECTANGLES_NUMBER; ++i)
	{
		Rectangle rectangle;
		rectangle.width = Common::Frand(RECTANGLES_SIZE_FROM, RECTANGLES_SIZE_TO);
		rectangle.height = Common::Frand(RECTANGLES_SIZE_FROM, RECTANGLES_SIZE_TO);
		rectangle.upperLeft.x = Common::Frand(0.0, RECTANGLES_SVG_WIDTH - rectangle.width);
		rectangle.upperLeft.y = Common::Frand(BORDER_ON_BOTTOM_AND_UPPER_PART + rectangle.height, 
			RECTANGLES_SVG_HEIGHT - BORDER_ON_BOTTOM_AND_UPPER_PART);

		problem.rectangles.push_back(std::move(rectangle));
	}

	// two explicit rectangles on borders
	static const double EXPLICIT_RECT_SIZE = 30.0;
	Rectangle leftRect{ { 0.0, RECTANGLES_SVG_HEIGHT/2.0 }, EXPLICIT_RECT_SIZE, EXPLICIT_RECT_SIZE };
	Rectangle rightRect{ { RECTANGLES_SVG_WIDTH - EXPLICIT_RECT_SIZE, RECTANGLES_SVG_HEIGHT / 2.0 }, EXPLICIT_RECT_SIZE, EXPLICIT_RECT_SIZE };
	problem.rectangles.push_back(leftRect);
	problem.rectangles.push_back(rightRect);

	return problem;
}

void VisualizeProblem(const char * name, const Problem & problem, Path * solution = nullptr)
{
	svg::Document doc(name,
		svg::Layout({ problem.width, problem.height }, svg::Layout::BottomLeft));

	for (const auto & rectangle : problem.rectangles)
	{
		doc << svg::Rectangle(svg::Point(rectangle.upperLeft.x, rectangle.upperLeft.y), 
			rectangle.width, rectangle.height, svg::Fill(svg::Color(100, 200, 120)));
	}

	if (solution)
	{
		std::vector<svg::Point> points(solution->size());
		for (size_t i = 0; i < solution->size(); ++i)
			points[i] = { solution->at(i).x, solution->at(i).y };

		doc << svg::Polyline(points, svg::Fill(), svg::Stroke(2.0, svg::Color(200, 100, 120)));
	}

	doc << svg::Rectangle({ 0.0, problem.height }, problem.width, problem.height,
		svg::Fill(), svg::Stroke(1.0, svg::Color::Silver));

	doc << svg::Circle(START_POINT, 6.0, svg::Color::Silver);

	doc.save();
}

Path ConvertToPath(const Problem & problem, const std::vector<double> & chromosome)
{
	Path ret;
	Common::Point previous, next;

	previous.x = chromosome[0];
	previous.y = chromosome[1];

	ret.push_back(START_POINT);

	previous.x += ret.back().x;
	previous.y += ret.back().y;

	if (ValidateLineSegment(problem, START_POINT, previous))
	{
		ret.push_back(previous);
		for (size_t i = 2; i < chromosome.size(); i = i + 2)
		{
			next.x = chromosome[i];
			next.y = chromosome[i+1];

			next.x += previous.x;
			next.y += previous.y;

			if (!ValidateLineSegment(problem, previous, next))
				break;

			previous = next;
			ret.push_back(previous);
		}
	}

	return ret;
}

double GetFitnessOfPath(const Path & path)
{
	return path.back().y;
	// penalize long path
	//return path.back().y - (double)path.size();
}

struct EvaluateRectangles
{
	BinaryGA::EvaluationResult operator()(uint32_t generation, const std::vector<double> & chromosome)
	{
		if (generation != currentGeneration)
		{
			std::cout << "\rGeneration " << std::fixed << generation;
			std::cout << " max height " << std::fixed << std::setprecision(2) << maxPath.back().y;
			//closestNumber = 0.0;
			currentGeneration = generation;
			numberOfGenerationsWithCurrentSolution++;
		}
		auto path = ConvertToPath(problem, chromosome);
		if (maxPath.empty() || GetFitnessOfPath(path) > GetFitnessOfPath(maxPath))
		{
			maxPath = path;
			numberOfGenerationsWithCurrentSolution = 0;

			return BinaryGA::EvaluationResult::ContinueProcessing;
		}

		if (numberOfGenerationsWithCurrentSolution >= STOP_AFTER_NUM_GENERATIONS_WITHOUT_CHANGE)
			return BinaryGA::EvaluationResult::ObjectiveReached;

		if (maxPath.back().y >= RECTANGLES_SVG_HEIGHT)
			return BinaryGA::EvaluationResult::ObjectiveReached;

		return BinaryGA::EvaluationResult::ContinueProcessing;
	}
	const Problem & problem;
	uint32_t currentGeneration = 0;
	uint32_t numberOfGenerationsWithCurrentSolution = 0;
	Path maxPath;
};


void TestRectangles()
{
	std::cout << "Rectangles" << std::endl;

	auto problem = GenerateProblem();
	VisualizeProblem("data\\RectanglesProblem.svg", problem);

	BinaryGA::Definition<double> definition;

	definition.parentSelection = BinaryGA::ParentSelectionType::Ranked;
	definition.mutation = BinaryGA::MutationType::Custom;
	definition.crossover = BinaryGA::CrossoverType::OnePoint;

	definition.populationSize = POPULATION_SIZE;
	definition.mutationProbability = MUTATION_PROBABILITY;
	definition.crossoverFactor = CROSSOVER_FACTOR;
	definition.maxNumberOfGenerations = MAX_NUMBER_OF_GENERATIONS;
	definition.numberOfGenes = MAXIMUM_PATH_SIZE * 2;

	definition.initializationCustomCallback = [](size_t index) -> std::vector<double>
	{
		std::vector<double> ret;
		double angle = (Common::PI / (double)POPULATION_SIZE) * (double)index;
		Common::Point point{ (double)STEP_MAX_LENGTH, 0.0 };
		point = Common::Rotate(point, angle);

		for (size_t i = 1; i < MAXIMUM_PATH_SIZE; ++i)
		{
			ret.push_back(point.x);
			ret.push_back(point.y);
		}

		return ret;
	};

	definition.computeFitness = [&problem](const std::vector<double> & chromosome) -> double
	{
		auto path = ConvertToPath(problem, chromosome);
		return GetFitnessOfPath(path);
	};

	definition.mutationCustomCallback = [](const double & value, size_t index) -> double
	{
		return Common::Frand(-STEP_MAX_LENGTH, STEP_MAX_LENGTH);
	};

	EvaluateRectangles evaluate{ problem };
	definition.evaluate = std::ref(evaluate);

	auto startTime = std::chrono::high_resolution_clock::now();
	auto solution = BinaryGA::Solve(definition);
	std::chrono::duration<double, std::milli> solveDuration = std::chrono::high_resolution_clock::now() - startTime;

	std::cout << std::endl << "Generation " << evaluate.currentGeneration << " (" << solveDuration.count() << "ms)" << std::endl;
	std::cout << "Best found solution: " << std::endl;
	std::cout << "height: " << std::fixed << std::setprecision(2) << evaluate.maxPath.back().y << std::endl;
	std::cout << std::endl;

	VisualizeProblem("data\\RectanglesSolution.svg", problem, &evaluate.maxPath);
}
