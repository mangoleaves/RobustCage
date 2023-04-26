#include "ColorMap.h"

namespace Cage
{
namespace SimpleUtils
{
// Default set maxValue to 1.0, minValue to 0.0
ColorMap::ColorMap()
{
	maxValue = 1.0;
	minValue = 0.0;
}

ColorMap::ColorMap(double maxV, double minV)
{
	maxValue = maxV;
	minValue = minV;
}

double ColorMap::GetMaxValue()
{
	return maxValue;
}

double ColorMap::GetMinValue()
{
	return minValue;
}

void ColorMap::SetMaxMinValue(double maxV, double minV)
{
	maxValue = maxV;
	minValue = minV;
}

Vec3d ColorMap::MapToColor(double value)
{
	double mapValue = (value - minValue) / (maxValue - minValue);
	if (mapValue < 0)
	{
		mapValue = 0;
	}
	else if (mapValue > 1)
	{
		mapValue = 1;
	}

	if (mapValue < 0.25)
	{
		return Vec3d(0, mapValue * 4, 1);
	}
	else if (mapValue < 0.5)
	{
		return Vec3d(0, 1, 1 - 4 * (mapValue - 0.25));
	}
	else if (mapValue < 0.75)
	{
		return Vec3d(4 * (mapValue - 0.5), 1, 0);
	}
	else
	{
		return Vec3d(1, 1 - 4 * (mapValue - 0.25), 0);
	}
}
}// namespace SimpleUtils
}// namespace Cage
