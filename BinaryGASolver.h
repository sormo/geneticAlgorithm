#include <vector>
#include <functional>

namespace BinaryGA
{
	enum class MutationType
	{
		None,
		Toggle,
		Swap,
		Custom
	};
	enum class CrossoverType
	{
		None,
		OnePoint,
		Ordered
	};
	enum class ParentSelectionType
	{
		RouletteWheel,
		Ranked
	};

	template<typename T>
	using ComputeFitnessCallback = std::function<double(const std::vector<T>&)>;
	enum EvaluationResult
	{
		ObjectiveReached,
		ContinueProcessing
	};
	template<typename T>
	using EvaluateCurrentState = std::function<EvaluationResult(uint32_t, const std::vector<T>&)>;

	// custom mutation callback receives value and index of gene in chromosome
	// return mutated value of gene
	template<typename T>
	using CustomMutation = std::function<T(const T&, size_t)>;

	template<typename T>
	using CustomInitialization = std::function<std::vector<T>(size_t)>;

	template<typename T>
	struct Definition
	{
		uint32_t numberOfGenes;
		uint32_t populationSize;
		uint32_t maxNumberOfGenerations;

		ComputeFitnessCallback<T> computeFitness;
		EvaluateCurrentState<T> evaluate;

		MutationType mutation;
		CrossoverType crossover;
		ParentSelectionType parentSelection;
		
		double mutationProbability;
		double crossoverFactor;

		CustomMutation<T> mutationCustomCallback;
		CustomInitialization<T> initializationCustomCallback;
	};

	template<typename T>
	std::vector<T> Solve(const Definition<T> & definition);

	// explicit template instantiation
	template std::vector<bool> Solve(const Definition<bool> & definition);
	template std::vector<uint8_t> Solve(const Definition<uint8_t> & definition);
	template std::vector<double> Solve(const Definition<double> & definition);
}
