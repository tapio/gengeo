#pragma once
#include <vector>
#include <cmath>
#include <functional>

namespace gengeo {

typedef unsigned int uint;

constexpr inline float abs(float a) { return a < 0 ? -a : a; }
constexpr inline float min(float a, float b) { return a < b ? a : b; }
constexpr inline float max(float a, float b) { return a > b ? a : b; }
constexpr inline float clamp(float x, float lo, float hi) { return x < lo ? lo : (x > hi ? hi : x); }
constexpr inline float saturate(float x) { return clamp(x, 0.f, 1.f); }
constexpr inline float lerp(float a, float b, float t) { return a + (b - a) * t; }
constexpr inline float sign(float x) { return x > 0.f ? 1.f : (x < 0.f ? -1.f : 0.f); }

struct vec2
{
	union {
		struct { float x, y; };
		struct { float r, g; };
		struct { float s, t; };
		float v[2];
	};

	constexpr vec2(float v = 0.f): x(v), y(v) {}
	constexpr vec2(float x_, float y_): x(x_), y(y_) {}

	constexpr vec2 operator+(const vec2& rhs) const { return { x + rhs.x, y + rhs.y }; }
	constexpr vec2 operator-(const vec2& rhs) const { return { x - rhs.x, y - rhs.y }; }
	constexpr vec2 operator*(const vec2& rhs) const { return { x * rhs.x, y * rhs.y }; }
	constexpr vec2 operator/(const vec2& rhs) const { return { x / rhs.x, y / rhs.y }; }
	constexpr vec2 operator*(float rhs) const { return { x * rhs, y * rhs }; }
	constexpr vec2 operator/(float rhs) const { return { x / rhs, y / rhs }; }

	constexpr vec2& operator+=(const vec2& rhs) { x += rhs.x; y += rhs.y; return *this; }
	constexpr vec2& operator-=(const vec2& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
	constexpr vec2& operator*=(const vec2& rhs) { x *= rhs.x; y *= rhs.y; return *this; }
	constexpr vec2& operator/=(const vec2& rhs) { x /= rhs.x; y /= rhs.y; return *this; }
	constexpr vec2& operator*=(float rhs) { x *= rhs; y *= rhs; return *this; }
	constexpr vec2& operator/=(float rhs) { x /= rhs; y /= rhs; return *this; }

	constexpr bool operator==(const vec2& rhs) const { return x == rhs.x && y == rhs.y; }
	constexpr bool operator!=(const vec2& rhs) const { return !(*this == rhs); }

	constexpr vec2& operator-() { x = -x; y = -y; return *this; }
	constexpr float operator[](int index) { return v[index]; }
};

struct vec3
{
	union {
		struct { float x, y, z; };
		struct { float r, g, b; };
		struct { float s, t, p; };
		float v[3];
	};

	constexpr vec3(float v = 0.f): x(v), y(v), z(v) {}
	constexpr vec3(float x_, float y_, float z_): x(x_), y(y_), z(z_) {}

	constexpr vec3 operator+(const vec3& rhs) const { return { x + rhs.x, y + rhs.y, z + rhs.z }; }
	constexpr vec3 operator-(const vec3& rhs) const { return { x - rhs.x, y - rhs.y, z - rhs.z }; }
	constexpr vec3 operator*(const vec3& rhs) const { return { x * rhs.x, y * rhs.y, z * rhs.z }; }
	constexpr vec3 operator/(const vec3& rhs) const { return { x / rhs.x, y / rhs.y, z / rhs.z }; }
	constexpr vec3 operator*(float rhs) const { return { x * rhs, y * rhs, z * rhs }; }
	constexpr vec3 operator/(float rhs) const { return { x / rhs, y / rhs, z / rhs }; }

	constexpr vec3& operator+=(const vec3& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
	constexpr vec3& operator-=(const vec3& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
	constexpr vec3& operator*=(const vec3& rhs) { x *= rhs.x; y *= rhs.y; z *= rhs.z; return *this; }
	constexpr vec3& operator/=(const vec3& rhs) { x /= rhs.x; y /= rhs.y; z /= rhs.z; return *this; }
	constexpr vec3& operator*=(float rhs) { x *= rhs; y *= rhs; z *= rhs; return *this; }
	constexpr vec3& operator/=(float rhs) { x /= rhs; y /= rhs; z /= rhs; return *this; }

	constexpr bool operator==(const vec3& rhs) const { return x == rhs.x && y == rhs.y && z == rhs.z; }
	constexpr bool operator!=(const vec3& rhs) const { return !(*this == rhs); }

	constexpr vec3& operator-() { x = -x; y = -y; z = -z; return *this; }
	constexpr float operator[](int index) { return v[index]; }
};

inline constexpr float dot(const vec2& lhs, const vec2& rhs) { return lhs.x * rhs.x + lhs.y * rhs.y; }
inline constexpr float length2(const vec2& v) { return dot(v, v); }
inline           float length(const vec2& v) { return std::sqrt(length2(v)); }
inline constexpr float distance2(const vec2& lhs, const vec2& rhs) { return length2(rhs - lhs); }
inline           float distance(const vec2& lhs, const vec2& rhs) { return length(rhs - lhs); }
inline           vec2 normalize(const vec2& v) { float l = 1.f / length(v); return v * l; }
inline constexpr vec2 abs(const vec2& v) { return { abs(v.x), abs(v.y) }; }
inline constexpr vec2 min(const vec2& lhs, const vec2& rhs) { return { min(lhs.x, rhs.x), min(lhs.y, rhs.y) }; }
inline constexpr vec2 max(const vec2& lhs, const vec2& rhs) { return { max(lhs.x, rhs.x), max(lhs.y, rhs.y) }; }
inline constexpr float min(const vec2& v) { return min(v.x, v.y); }
inline constexpr float max(const vec2& v) { return max(v.x, v.y); }
inline constexpr vec2 clamp(const vec2& v, const vec2& lo, const vec2& hi) { return { clamp(v.x, lo.x, hi.x), clamp(v.y, lo.y, hi.y) }; }
inline constexpr vec2 lerp(const vec2& a, const vec2& b, float t) { return a + (b - a) * t; }

inline constexpr vec3 cross(const vec3& lhs, const vec3& rhs) { return { lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x }; }
inline constexpr float dot(const vec3& lhs, const vec3& rhs) { return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z; }
inline constexpr float length2(const vec3& v) { return dot(v, v); }
inline           float length(const vec3& v) { return std::sqrt(length2(v)); }
inline constexpr float distance2(const vec3& lhs, const vec3& rhs) { return length2(rhs - lhs); }
inline           float distance(const vec3& lhs, const vec3& rhs) { return length(rhs - lhs); }
inline           vec3 normalize(const vec3& v) { float l = 1.f / length(v); return v * l; }
inline constexpr vec3 abs(const vec3& v) { return { abs(v.x), abs(v.y), abs(v.z) }; }
inline constexpr vec3 min(const vec3& lhs, const vec3& rhs) { return { min(lhs.x, rhs.x), min(lhs.y, rhs.y), min(lhs.z, rhs.z) }; }
inline constexpr vec3 max(const vec3& lhs, const vec3& rhs) { return { max(lhs.x, rhs.x), max(lhs.y, rhs.y), max(lhs.z, rhs.z) }; }
inline constexpr float min(const vec3& v) { return min(v.x, min(v.y, v.z)); }
inline constexpr float max(const vec3& v) { return max(v.x, max(v.y, v.z)); }
inline constexpr vec3 clamp(const vec3& v, const vec3& lo, const vec3& hi) { return { clamp(v.x, lo.x, hi.x), clamp(v.y, lo.y, hi.y), clamp(v.z, lo.z, hi.z) }; }
inline constexpr vec3 lerp(const vec3& a, const vec3& b, float t) { return a + (b - a) * t; }
inline           vec3 triangleNormal(const vec3& a, const vec3& b, const vec3& c) { return normalize(cross(b - a, c - a)); }

struct Bounds {
	vec3 min;
	vec3 max;

	constexpr vec3 center() const { return (min + max) * 0.5f; }
	constexpr vec3 size() const { return max - min; }
	constexpr vec3 extents() const { return (max - min) * 0.5f; }
	constexpr Bounds& expand(const vec3& amount ) { min -= amount; max += amount; return *this; }
	constexpr bool contains(const vec3& p) const {
		return min.x <= p.x && min.y <= p.y && min.z <= p.y && max.x >= p.x && max.y >= p.y && max.z >= p.z;
	}
};


struct TriFace
{
	uint a, b, c;
	constexpr TriFace(uint a_, uint b_, uint c_): a(a_), b(b_), c(c_) {}
	constexpr void offset(int offset) { a += offset; b += offset; c += offset; }
};

struct QuadFace
{
	uint a, b, c, d;
	constexpr QuadFace(uint a_, uint b_, uint c_, uint d_): a(a_), b(b_), c(c_), d(d_) {}
	constexpr void offset(int offset) { a += offset; b += offset; c += offset; d += offset; }
};

struct Geometry
{
	std::vector<TriFace> triangles;
	std::vector<QuadFace> quads;
	std::vector<vec3> positions;
	std::vector<vec2> texcoords;
	std::vector<vec3> normals;
};

void bounds(const Geometry& geo, vec3& lower, vec3& upper);
Bounds bounds(const Geometry& geo);
float boundingSphere(const Geometry& geo);

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

// Unit cube
Geometry cube();
// Sized box
inline Geometry box(vec3 size) { auto c = cube(); return scale(c, size); }
// Unit sphere
Geometry sphere(int subdivs = 4);
// Unit circle
Geometry circley(int points = 16);
// "Planes" are unit squares, with normal pointing to positive
Geometry planexy();
Geometry planexz();
Geometry planeyz();
inline Geometry planex() { return planeyz(); }
inline Geometry planey() { return planexz(); }
inline Geometry planez() { return planexy(); }

Geometry heightmap(float* heights, int w, int h, vec3 scale = vec3(1));

Geometry& subdivide(Geometry& geo, int amount = 1);
Geometry& weld(Geometry& geo, float maxNormalAngleDeg = 30.f, double tolerance = 0.0001);
Geometry& triangulate(Geometry& geo);
Geometry& displaceAlongNormals(Geometry& geo, float minAmount, float maxAmount, int seed = 0);
Geometry& displacementMap(Geometry& geo, float* displacementTexture, int w, int h, float height = 1.0f);
Geometry& calculateNormals(Geometry& geo);
Geometry& normalizeNormals(Geometry& geo);
Geometry& flipNormals(Geometry& geo);

bool writeObj(const Geometry& geo, const char* filename = "out.obj");


namespace sdf {

	inline float sphere(vec3 p, float r) {
		return length(p) - r;
	}
	inline float box(vec3 p, vec3 b) {
		vec3 d = abs(p) - b;
		return length(max(d, vec3(0.0f))) + min(max(d), 0.0f);
	}
	inline constexpr float plane(vec3 p, vec3 n, float d) {
		return dot(p, n) + d; // n must be normalized
	}
	inline float capsule(vec3 p, vec3 a, vec3 b, float r) {
		vec3 pa = p - a, ba = b - a;
		float h = saturate(dot(pa, ba) / dot(ba, ba));
		return length(pa - ba * h) - r;
	}
	inline float capsulex(vec3 p, float h, float r) {
		p.x -= clamp(p.x, 0.0f, h);
		return length(p) - r;
	}
	inline float capsuley(vec3 p, float h, float r) {
		p.y -= clamp(p.y, 0.0f, h);
		return length(p) - r;
	}
	inline float capsulez(vec3 p, float h, float r) {
		p.z -= clamp(p.z, 0.0f, h);
		return length(p) - r;
	}
	inline float cylinder(vec3 p, vec3 a, vec3 b, float r) {
		vec3  ba = b - a;
		vec3  pa = p - a;
		float baba = dot(ba, ba);
		float paba = dot(pa, ba);
		float x = length(pa * baba - ba * paba) - r * baba;
		float y = abs(paba - baba * 0.5f) - baba * 0.5f;
		float x2 = x * x;
		float y2 = y * y * baba;
		float d = (max(x, y) < 0.0f) ? -min(x2, y2): (((x > 0.0f) ? x2 : 0.0f) + ((y > 0.0f) ? y2 : 0.0f));
		return sign(d) * sqrt(abs(d)) / baba;
	}
	inline float cylinderx(vec3 p, float h, float r) {
		vec2 d = abs(vec2(length(vec2(p.y, p.z)), p.x)) - vec2(h, r);
		return min(max(d.x, d.y), 0.0f) + length(max(d, 0.0f));
	}
	inline float cylindery(vec3 p, float h, float r) {
		vec2 d = abs(vec2(length(vec2(p.x, p.z)), p.y)) - vec2(h, r);
		return min(max(d.x, d.y), 0.0f) + length(max(d, 0.0f));
	}
	inline float cylinderz(vec3 p, float h, float r) {
		vec2 d = abs(vec2(length(vec2(p.x, p.y)), p.z)) - vec2(h, r);
		return min(max(d.x, d.y), 0.0f) + length(max(d, 0.0f));
	}


	inline constexpr float opRound(float sdf, float r) { return sdf - r; }
	inline constexpr float opUnion(float a, float b) { return min(a, b); }
	inline constexpr float opSubtract(float a, float b) { return max(a, -b); }
	inline constexpr float opIntersect(float a, float b) { return max(a, b); }
	inline constexpr float opSmoothUnion(float a, float b, float k) {
		const float h = saturate(0.5f + 0.5f * (b - a) / k);
		return lerp(b, a, h) - k * h * (1.0f - h);
	}
	inline constexpr float opSmoothSubtraction(float a, float b, float k) {
		const float h = saturate(0.5f - 0.5f * (b + a) / k);
		return lerp(-b, a, h) + k * h * (1.0f - h);
	}
	inline constexpr float opSmoothIntersection(float a, float b, float k) {
		const float h = saturate(0.5f - 0.5f * (b - a) / k);
		return lerp(b, a, h) + k * h * (1.0f - h);
	}

} // namespace sdf

typedef std::function<float(vec3)> SDFFunction;
Geometry polygonize(vec3 start, vec3 end, vec3 step, SDFFunction sdf);



struct VoxelGrid
{
	typedef char cell_t;

	VoxelGrid(int xsize, int ysize, int zsize): xsize(xsize), ysize(ysize), zsize(zsize) {
		cells.resize(xsize * ysize * zsize);
	}

	cell_t get(int x, int y, int z) const { return cells[cellIndex(x, y, z)]; }
	cell_t get(vec3 pos) const { return get(pos.x, pos.y, pos.z); }
	cell_t getSafe(int x, int y, int z, cell_t def = cell_t()) const {
		return (x >= 0 && x < xsize && y >= 0 && y < ysize && z >= 0 && z < zsize) ? get(x, y, z) : def; }
	cell_t getSafe(vec3 pos, cell_t def = cell_t()) const { return getSafe(pos.x, pos.y, pos.z, def); }
	
	void set(int x, int y, int z, const cell_t& cell) { cells[cellIndex(x, y, z)] = cell; }
	void set(vec3 pos, const cell_t& cell) { set(pos.x, pos.y, pos.z, cell); }
	
	int cellIndex(int x, int y, int z) const { return (z * xsize * ysize) + (y * xsize) + x; }

	template<typename VisitFunc>
	void visit(VisitFunc visitFunc) const {
		for (int k = 0; k < zsize; ++k) {
			for (int j = 0; j < ysize; ++j) {
				for (int i = 0; i < xsize; ++i) {
					visitFunc(i, j, k, get(i, j, k));
				}
			}
		}
	}

	template<typename VisitFunc>
	void visit(VisitFunc visitFunc) {
		for (int k = 0; k < zsize; ++k) {
			for (int j = 0; j < ysize; ++j) {
				for (int i = 0; i < xsize; ++i) {
					set(i, j, k, visitFunc(i, j, k));
				}
			}
		}
	}

	template<typename VisitFunc>
	void visit(vec3 start, vec3 end, VisitFunc visitFunc) {
		for (int k = start.z; k < end.z; ++k) {
			for (int j = start.y; j < end.y; ++j) {
				for (int i = start.x; i < end.x; ++i) {
					set(i, j, k, visitFunc(i, j, k));
				}
			}
		}
	}

	int xsize = 0;
	int ysize = 0;
	int zsize = 0;

	std::vector<cell_t> cells;
};

VoxelGrid& box(VoxelGrid& voxelGrid, vec3 start, vec3 end, VoxelGrid::cell_t cell = 1);
Geometry polygonize(const VoxelGrid& voxelGrid);

Geometry polygonizeSmooth(const VoxelGrid& voxelGrid, float isolevel = 0);

inline void hash_combine(std::size_t& seed, std::size_t v) {
	seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

} // namespace gengeo

namespace std {
	template<> struct hash<gengeo::vec2> {
		size_t operator()(const gengeo::vec2& v) const {
			size_t seed = 34598591;
			std::hash<float> hasher;
			gengeo::hash_combine(seed, hasher(v.x));
			gengeo::hash_combine(seed, hasher(v.y));
			return seed;
		}
	};
	template<> struct hash<gengeo::vec3> {
		size_t operator()(const gengeo::vec3& v) const {
			size_t seed = 81527137;
			std::hash<float> hasher;
			gengeo::hash_combine(seed, hasher(v.x));
			gengeo::hash_combine(seed, hasher(v.y));
			gengeo::hash_combine(seed, hasher(v.z));
			return seed;
		}
	};
}