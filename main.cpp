#include <cstdint>
#include <ctime>
#include <iostream>

void TestKnapsack();
void TestDigitsAndOperators();
void TestTravelingSalesman();
void TestCircles();
void TestRectangles();

int main()
{
	srand((int)time(NULL));

	TestKnapsack();
	TestDigitsAndOperators();
	TestTravelingSalesman();
	TestCircles();
	TestRectangles();

	return 0;
}
