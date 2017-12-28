#include "Common.h"
#include <cstdlib>
#include <vector>

namespace Common
{
	double Frand(double min, double max)
	{
		double f = (double)rand() / RAND_MAX;
		return min + f * (max - min);
	}

	double Distance(const Point & first, const Point & second)
	{
		return sqrt(pow(first.x - second.x, 2) + pow(first.y - second.y, 2));
	}

	uint64_t ConvertBitsAtIndexToIntegralNumber(const std::vector<bool> & chromosome, size_t index)
	{
		uint64_t ret = 0;
		for (size_t i = 0; i < NUMBER_OF_BITS_IN_DOUBLE; ++i)
		{
			if (chromosome[index + i])
				ret |= (1ULL << i);
		}
		return ret;
	}

	double ConvertIntegralNumberToDoubleInRange(uint64_t number, double from, double to)
	{
		return ((double)number / (double)std::numeric_limits<uint64_t>::max())*(double)to + from;
	}

	Point Rotate(const Point & v, double angle)
	{
		double s = std::sin(angle);
		double c = std::cos(angle);

		Point t;

		t.x = v.x * c - v.y * s;
		t.y = v.x * s + v.y * c;

		return t;
	}
}
