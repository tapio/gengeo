#pragma once
#include <vector>
#include <cmath>

namespace gengeo {

typedef unsigned int uint;

using std::abs;
using std::min;
using std::max;

struct vec2
{
	union {
		struct { float x, y; };
		struct { float r, g; };
		struct { float s, t; };
	};

	vec2(float v = 0.f): x(v), y(v) {}
	vec2(float x_, float y_): x(x_), y(y_) {}

	vec2 operator+(const vec2& rhs) const { return vec2(x + rhs.x, y + rhs.y); }
	vec2 operator-(const vec2& rhs) const { return vec2(x - rhs.x, y - rhs.y); }
	vec2 operator*(const vec2& rhs) const { return vec2(x * rhs.x, y * rhs.y); }
	vec2 operator/(const vec2& rhs) const { return vec2(x / rhs.x, y / rhs.y); }
	vec2 operator*(float rhs) const { return vec2(x * rhs, y * rhs); }
	vec2 operator/(float rhs) const { return vec2(x / rhs, y / rhs); }

	vec2& operator+=(const vec2& rhs) { x += rhs.x; y += rhs.y; return *this; }
	vec2& operator-=(const vec2& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
	vec2& operator*=(const vec2& rhs) { x *= rhs.x; y *= rhs.y; return *this; }
	vec2& operator/=(const vec2& rhs) { x /= rhs.x; y /= rhs.y; return *this; }
	vec2& operator*=(float rhs) { x *= rhs; y *= rhs; return *this; }
	vec2& operator/=(float rhs) { x /= rhs; y /= rhs; return *this; }

	vec2& operator-() { x = -x; y = -y; return *this; }
};

struct vec3
{
	union {
		struct { float x, y, z; };
		struct { float r, g, b; };
		struct { float s, t, p; };
	};

	vec3(float v = 0.f): x(v), y(v), z(v) {}
	vec3(float x_, float y_, float z_): x(x_), y(y_), z(z_) {}

	vec3 operator+(const vec3& rhs) const { return vec3(x + rhs.x, y + rhs.y, z + rhs.z); }
	vec3 operator-(const vec3& rhs) const { return vec3(x - rhs.x, y - rhs.y, z - rhs.z); }
	vec3 operator*(const vec3& rhs) const { return vec3(x * rhs.x, y * rhs.y, z * rhs.z); }
	vec3 operator/(const vec3& rhs) const { return vec3(x / rhs.x, y / rhs.y, z / rhs.z); }
	vec3 operator*(float rhs) const { return vec3(x * rhs, y * rhs, z * rhs); }
	vec3 operator/(float rhs) const { return vec3(x / rhs, y / rhs, z / rhs); }

	vec3& operator+=(const vec3& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
	vec3& operator-=(const vec3& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
	vec3& operator*=(const vec3& rhs) { x *= rhs.x; y *= rhs.y; z *= rhs.z; return *this; }
	vec3& operator/=(const vec3& rhs) { x /= rhs.x; y /= rhs.y; z /= rhs.z; return *this; }
	vec3& operator*=(float rhs) { x *= rhs; y *= rhs; z *= rhs; return *this; }
	vec3& operator/=(float rhs) { x /= rhs; y /= rhs; z /= rhs; return *this; }

	vec3& operator-() { x = -x; y = -y; z = -z; return *this; }
};

inline float dot(const vec2& lhs, const vec2& rhs) { return lhs.x * rhs.x + lhs.y * rhs.y; }
inline float length2(const vec2& v) { return dot(v, v); }
inline float length(const vec2& v) { return std::sqrt(length2(v)); }
inline float distance2(const vec2& lhs, const vec2& rhs) { return length2(rhs - lhs); }
inline float distance(const vec2& lhs, const vec2& rhs) { return length(rhs - lhs); }
inline vec2 normalize(const vec2& v) { float l = 1.f / length(v); return v * l; }
inline vec2 abs(const vec2& v) { return vec2(std::abs(v.x), std::abs(v.y)); }
inline vec2 min(const vec2& lhs, const vec2& rhs) { return vec2(std::min(lhs.x, rhs.x), std::min(lhs.y, rhs.y)); }
inline vec2 max(const vec2& lhs, const vec2& rhs) { return vec2(std::max(lhs.x, rhs.x), std::max(lhs.y, rhs.y)); }

inline vec3 cross(const vec3& lhs, const vec3& rhs) { return vec3(lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x); }
inline float dot(const vec3& lhs, const vec3& rhs) { return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z; }
inline float length2(const vec3& v) { return dot(v, v); }
inline float length(const vec3& v) { return std::sqrt(length2(v)); }
inline float distance2(const vec3& lhs, const vec3& rhs) { return length2(rhs - lhs); }
inline float distance(const vec3& lhs, const vec3& rhs) { return length(rhs - lhs); }
inline vec3 normalize(const vec3& v) { float l = 1.f / length(v); return v * l; }
inline vec3 abs(const vec3& v) { return vec3(std::abs(v.x), std::abs(v.y), std::abs(v.z)); }
inline vec3 min(const vec3& lhs, const vec3& rhs) { return vec3(std::min(lhs.x, rhs.x), std::min(lhs.y, rhs.y), std::min(lhs.z, rhs.z)); }
inline vec3 max(const vec3& lhs, const vec3& rhs) { return vec3(std::max(lhs.x, rhs.x), std::max(lhs.y, rhs.y), std::max(lhs.z, rhs.z)); }


struct TriFace
{
	uint a, b, c;
	TriFace(uint a_, uint b_, uint c_): a(a_), b(b_), c(c_) {}
	void offset(int offset) { a += offset; b += offset; c += offset; }
};

struct QuadFace
{
	uint a, b, c, d;
	QuadFace(uint a_, uint b_, uint c_, uint d_): a(a_), b(b_), c(c_), d(d_) {}
	void offset(int offset) { a += offset; b += offset; c += offset; d += offset; }
};

struct Geometry
{
	std::vector<TriFace> triangles;
	std::vector<QuadFace> quads;
	std::vector<vec3> positions;
	std::vector<vec2> texcoords;
	std::vector<vec3> normals;
};

// Unit cube
Geometry cube();
// Unit sphere
Geometry sphere(int subdivs = 4);
// "Planes" are unit squares, with normal pointing to positive
Geometry planexy();
Geometry planexz();
Geometry planeyz();

vec3 extents(const Geometry& geo);
void bounds(const Geometry& geo, vec3& lower, vec3& upper);

Geometry& scale(Geometry& geo, vec3 s);
Geometry& translate(Geometry& geo, vec3 v);
Geometry& center(Geometry& geo);
Geometry& concat(Geometry& geo, const Geometry& other);
inline Geometry& scalex(Geometry& geo, float s) { return scale(geo, {s, 1, 1}); }
inline Geometry& scaley(Geometry& geo, float s) { return scale(geo, {1, s, 1}); }
inline Geometry& scalez(Geometry& geo, float s) { return scale(geo, {1, 1, s}); }
inline Geometry& mirrorx(Geometry& geo) { return scalex(geo, -1); }
inline Geometry& mirrory(Geometry& geo) { return scaley(geo, -1); }
inline Geometry& mirrorz(Geometry& geo) { return scalez(geo, -1); }
inline Geometry& translatex(Geometry& geo, float x) { return translate(geo, {x, 0, 0}); }
inline Geometry& translatey(Geometry& geo, float y) { return translate(geo, {0, y, 0}); }
inline Geometry& translatez(Geometry& geo, float z) { return translate(geo, {0, 0, z}); }

Geometry& subdivide(Geometry& geo, int amount = 1);
Geometry& triangulate(Geometry& geo);
Geometry& normalizeNormals(Geometry& geo);

bool writeObj(const Geometry& geo, const char* filename = "out.obj");


} // namespace gengeo
