#pragma once
#include <cstdint>
#include <vector>

#define NUMBER_OF_BITS_IN_DOUBLE (sizeof(double)*8)

namespace Common
{
	static const double PI = 3.141592654;

	struct Point
	{
		double x;
		double y;
	};

	double Frand(double min, double max);
	double Distance(const Point & first, const Point & second);
	Point Rotate(const Point & v, double angle);

	uint64_t ConvertBitsAtIndexToIntegralNumber(const std::vector<bool> & chromosome, size_t index);
	double ConvertIntegralNumberToDoubleInRange(uint64_t number, double from, double to);
}