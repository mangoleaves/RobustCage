#pragma once

#include <cmath>
#include <vector>
#include <array>

namespace Cage
{
namespace Geometry
{
class Vec2d
{
public:
	double _x;
	double _y;
public:
	Vec2d() {}
	Vec2d(const double& v0, const double& v1) :_x(v0), _y(v1){}
	Vec2d(const double* v) :_x(v[0]), _y(v[1]) {}
	Vec2d(const Vec2d& v) :_x(v._x), _y(v._y) {}

	Vec2d& operator=(const Vec2d& v)
	{
		_x = v._x; _y = v._y;
		return *this;
	}

	~Vec2d() {}

	inline double& x(){ return _x; }
	inline double& y(){ return _y; }
	inline double x()const { return _x; }
	inline double y()const  { return _y; }

	double* data() { return reinterpret_cast<double*>(this); }
	const double* data() const { return reinterpret_cast<const double*>(this); }

	double& operator[](size_t dim) { return reinterpret_cast<double*>(this)[dim]; }
	const double& operator[](size_t dim)const { return reinterpret_cast<const double*>(this)[dim]; }

	inline Vec2d operator-() const { return Vec2d(-_x, -_y); }
	inline Vec2d operator-(const Vec2d& rhs) const { return Vec2d(_x - rhs._x, _y - rhs._y); }
	inline Vec2d operator+(const Vec2d& rhs) const { return Vec2d(_x + rhs._x, _y + rhs._y); }
	inline Vec2d operator*(const double rhs) const { return Vec2d(_x * rhs, _y * rhs); }
	inline Vec2d operator/(const double rhs) const { return Vec2d(_x / rhs, _y / rhs); }

	inline bool operator==(const Vec2d& rhs) const { return _x == rhs._x && _y == rhs._y; }
	inline bool operator!=(const Vec2d& rhs) const { return !(*this == rhs); }
	inline bool less_on(size_t dim, const Vec2d& rhs)const { return operator[](dim) < rhs[dim]; }
	inline bool less_on(size_t dim, double rhs)const { return operator[](dim) < rhs; }

	inline double dot(const Vec2d& rhs) const { return _x * rhs._x + _y * rhs._y; }

	inline Vec2d normalized() const
	{
		double z = squaredNorm();
		return z > 0 ? (*this) / sqrt(z) : (*this);
	}

	inline double squaredNorm() const { return _x * _x + _y * _y; }
	inline double norm() const { return std::sqrt(squaredNorm()); }

	inline void minimize(const Vec2d& rhs)
	{
		if (rhs._x < _x)_x = rhs._x;
		if (rhs._y < _y)_y = rhs._y;
	}
	inline void maximize(const Vec2d& rhs)
	{
		if (rhs._x > _x)_x = rhs._x;
		if (rhs._y > _y)_y = rhs._y;
	}
};

inline Vec2d operator*(const double lhs, const Vec2d& rhs) { return rhs * lhs; }
}
// namespace Geometry
}// namespace Cage