#pragma once
#include "Geometry/Basic/Types.h"

namespace Cage
{
namespace SimpleUtils
{

using namespace Geometry;

class ColorMap
{
private:
	double maxValue, minValue;
public:
	ColorMap();
	ColorMap(double maxV, double minV);

	double GetMaxValue();
	double GetMinValue();
	void SetMaxMinValue(double maxV, double minV);

	Vec3d MapToColor(double value);
};
}// namespace SimpleUtils
}// namespace Cage

