#pragma once
#include <vector>
#include <cmath>
#include <functional>

namespace gengeo {

static constexpr float PI = 3.1415926535f;
static constexpr float TWOPI = 2 * PI;
//static constexpr float PI2 = PI * 0.5f;
//static constexpr float PI4 = PI * 0.25f;
static constexpr float DEG_TO_RAD = PI / 180.f;
static constexpr float RAD_TO_DEG = 180.f / PI;

typedef unsigned int uint;

template<typename T>
constexpr inline T abs(T a) { return a < 0 ? -a : a; }
template<typename T>
constexpr inline T min(T a, T b) { return a < b ? a : b; }
template<typename T>
constexpr inline T max(T a, T b) { return a > b ? a : b; }
template<typename T>
constexpr inline T clamp(T x, T lo, T hi) { return x < lo ? lo : (x > hi ? hi : x); }
constexpr inline float saturate(float x) { return clamp(x, 0.f, 1.f); }
constexpr inline float fract(float x) { return x - (long)x; }
constexpr inline int mod(int a, int b) { int c = a % b; return (c < 0) ? c + b : c; }
constexpr inline float fmod(float a, float b) { float c = a - (int)(a / b) * b; return (c < 0) ? c + b : c; }
constexpr inline float lerp(float a, float b, float t) { return a + (b - a) * t; }
constexpr inline float sign(float x) { return x > 0.f ? 1.f : (x < 0.f ? -1.f : 0.f); }
constexpr inline float step(float edge, float x) { return x < edge ? 0.f : 1.f; }
constexpr inline float smoothstep(float edge0, float edge1, float x) { float t = clamp((x - edge0) / (edge1 - edge0), 0.0f, 1.0f); return t * t * (3.0f - 2.0f * t); }
constexpr inline int roundToInt(float x) { return (int)(x + 0.5f - (x < 0)); }
constexpr inline float round(float x) { return (float)roundToInt(x); }
constexpr inline float mapRange(float value, float inMin, float inMax, float outMin, float outMax) { return (value - inMin) / (inMax - inMin) * (outMax - outMin) + outMin; }
constexpr inline float mapRangeClamp(float value, float inMin, float inMax, float outMin, float outMax) { return (clamp(value, inMin, inMax) - inMin) / (inMax - inMin) * (outMax - outMin) + outMin; }

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

	constexpr vec2 operator-() const { return { -x, -y }; }
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

	constexpr float& operator[](int index) { return v[index]; }
	constexpr float operator[](int index) const { return v[index]; }
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

	constexpr vec3 operator-() const { return { -x, -y, -z }; }
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

	constexpr float& operator[](int index) { return v[index]; }
	constexpr float operator[](int index) const { return v[index]; }

	constexpr vec2 xy() const { return { x, y }; }
};

struct vec4
{
	union {
		struct { float x, y, z, w; };
		struct { float r, g, b, a; };
		struct { float s, t, p, q; };
		// GCC doesn't like these unions because of constructors...
		//struct { vec3 xyz; float padding1; };
		//struct { vec2 xy; vec2 zw; };
		float v[4];
	};

	constexpr vec4(float v = 0.f) : x(v), y(v), z(v), w(v) {}
	constexpr vec4(float x_, float y_, float z_, float w_) : x(x_), y(y_), z(z_), w(w_) {}

	constexpr vec4 operator-() const { return { -x, -y, -z, -w }; }
	constexpr vec4 operator+(const vec4& rhs) const { return { x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w }; }
	constexpr vec4 operator-(const vec4& rhs) const { return { x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w }; }
	constexpr vec4 operator*(const vec4& rhs) const { return { x * rhs.x, y * rhs.y, z * rhs.z, w * rhs.w }; }
	constexpr vec4 operator/(const vec4& rhs) const { return { x / rhs.x, y / rhs.y, z / rhs.z, w / rhs.w }; }
	constexpr vec4 operator*(float rhs) const { return { x * rhs, y * rhs, z * rhs, w * rhs }; }
	constexpr vec4 operator/(float rhs) const { return { x / rhs, y / rhs, z / rhs, w / rhs }; }

	constexpr vec4& operator+=(const vec4& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; w += rhs.w; return *this; }
	constexpr vec4& operator-=(const vec4& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; w -= rhs.w; return *this; }
	constexpr vec4& operator*=(const vec4& rhs) { x *= rhs.x; y *= rhs.y; z *= rhs.z; w *= rhs.w; return *this; }
	constexpr vec4& operator/=(const vec4& rhs) { x /= rhs.x; y /= rhs.y; z /= rhs.z; w /= rhs.w; return *this; }
	constexpr vec4& operator*=(float rhs) { x *= rhs; y *= rhs; z *= rhs; w *= rhs; return *this; }
	constexpr vec4& operator/=(float rhs) { x /= rhs; y /= rhs; z /= rhs; w /= rhs; return *this; }

	constexpr bool operator==(const vec4& rhs) const { return x == rhs.x && y == rhs.y && z == rhs.z && w == rhs.w; }
	constexpr bool operator!=(const vec4& rhs) const { return !(*this == rhs); }

	constexpr float& operator[](int index) { return v[index]; }
	constexpr float operator[](int index) const { return v[index]; }

	constexpr vec2 xy() const { return { x, y }; }
	constexpr vec3 xyz() const { return { x, y, z }; }
};

constexpr vec2 operator*(float lhs, const vec2& rhs) { return rhs * lhs; }
constexpr vec3 operator*(float lhs, const vec3& rhs) { return rhs * lhs; }
constexpr vec4 operator*(float lhs, const vec4& rhs) { return rhs * lhs; }

inline constexpr float dot(const vec2& lhs, const vec2& rhs) { return lhs.x * rhs.x + lhs.y * rhs.y; }
inline constexpr float length2(const vec2& v) { return dot(v, v); }
inline           float length(const vec2& v) { return std::sqrt(length2(v)); }
inline constexpr float distance2(const vec2& lhs, const vec2& rhs) { return length2(rhs - lhs); }
inline           float distance(const vec2& lhs, const vec2& rhs) { return length(rhs - lhs); }
inline           vec2 normalize(const vec2& v) { float l = 1.f / length(v); return v * l; }
inline constexpr vec2 round(const vec2& v) { return { round(v.x), round(v.y) }; }
inline constexpr vec2 abs(const vec2& v) { return { abs(v.x), abs(v.y) }; }
inline constexpr vec2 min(const vec2& lhs, const vec2& rhs) { return { min(lhs.x, rhs.x), min(lhs.y, rhs.y) }; }
inline constexpr vec2 max(const vec2& lhs, const vec2& rhs) { return { max(lhs.x, rhs.x), max(lhs.y, rhs.y) }; }
inline constexpr float min(const vec2& v) { return min(v.x, v.y); }
inline constexpr float max(const vec2& v) { return max(v.x, v.y); }
inline constexpr vec2 clamp(const vec2& v, const vec2& lo, const vec2& hi) { return { clamp(v.x, lo.x, hi.x), clamp(v.y, lo.y, hi.y) }; }
inline constexpr vec2 lerp(const vec2& a, const vec2& b, float t) { return a + (b - a) * t; }
inline constexpr vec2 lerp(const vec2& a, const vec2& b, const vec2& t) { return { lerp(a.x, b.x, t.x), lerp(a.y, b.y, t.y) }; }

inline constexpr vec3 cross(const vec3& lhs, const vec3& rhs) { return { lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x }; }
inline constexpr float dot(const vec3& lhs, const vec3& rhs) { return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z; }
inline constexpr float length2(const vec3& v) { return dot(v, v); }
inline           float length(const vec3& v) { return std::sqrt(length2(v)); }
inline constexpr float distance2(const vec3& lhs, const vec3& rhs) { return length2(rhs - lhs); }
inline           float distance(const vec3& lhs, const vec3& rhs) { return length(rhs - lhs); }
inline           vec3 normalize(const vec3& v) { float l = 1.f / length(v); return v * l; }
inline           vec3 orthonormalize(const vec3& forward, vec3& up) { vec3 right = normalize(cross(up, forward)); up = normalize(cross(forward, right)); return right; }
inline constexpr vec3 round(const vec3& v) { return { round(v.x), round(v.y), round(v.z) }; }
inline constexpr vec3 abs(const vec3& v) { return { abs(v.x), abs(v.y), abs(v.z) }; }
inline constexpr vec3 min(const vec3& lhs, const vec3& rhs) { return { min(lhs.x, rhs.x), min(lhs.y, rhs.y), min(lhs.z, rhs.z) }; }
inline constexpr vec3 max(const vec3& lhs, const vec3& rhs) { return { max(lhs.x, rhs.x), max(lhs.y, rhs.y), max(lhs.z, rhs.z) }; }
inline constexpr float min(const vec3& v) { return min(v.x, min(v.y, v.z)); }
inline constexpr float max(const vec3& v) { return max(v.x, max(v.y, v.z)); }
inline constexpr vec3 clamp(const vec3& v, const vec3& lo, const vec3& hi) { return { clamp(v.x, lo.x, hi.x), clamp(v.y, lo.y, hi.y), clamp(v.z, lo.z, hi.z) }; }
inline           float uangle(const vec3& a, const vec3& b) { float d = dot(a, b); return d > 1 ? 0 : std::acos(d < -1 ? -1 : d); }
inline           float angle(const vec3& a, const vec3& b) { return uangle(normalize(a), normalize(b)); }
inline constexpr vec3 lerp(const vec3& a, const vec3& b, float t) { return a + (b - a) * t; }
inline constexpr vec3 lerp(const vec3& a, const vec3& b, const vec3& t) { return { lerp(a.x, b.x, t.x), lerp(a.y, b.y, t.y), lerp(a.z, b.z, t.z) }; }
inline           vec3 nlerp(const vec3& a, const vec3& b, float t) { return normalize(lerp(a, b, t)); }
inline           vec3 slerp(const vec3& a, const vec3& b, float t) { float th = uangle(a, b); return th == 0 ? a : a * (std::sin(th * (1 - t)) / std::sin(th)) + b * (std::sin(th * t) / std::sin(th)); }
inline           vec3 triangleNormal(const vec3& a, const vec3& b, const vec3& c) { return normalize(cross(b - a, c - a)); }
inline           vec3 quadNormal(const vec3& a, const vec3& b, const vec3& c, const vec3& d) { return normalize(cross(c - a, d - b)); }

inline constexpr float dot(const vec4& lhs, const vec4& rhs) { return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w; }
inline constexpr float length2(const vec4& v) { return dot(v, v); }
inline           float length(const vec4& v) { return std::sqrt(length2(v)); }
inline constexpr float distance2(const vec4& lhs, const vec4& rhs) { return length2(rhs - lhs); }
inline           float distance(const vec4& lhs, const vec4& rhs) { return length(rhs - lhs); }
inline           vec4 normalize(const vec4& v) { float l = 1.f / length(v); return v * l; }
inline constexpr vec4 round(const vec4& v) { return { round(v.x), round(v.y), round(v.z), round(v.w) }; }
inline constexpr vec4 abs(const vec4& v) { return { abs(v.x), abs(v.y), abs(v.z), abs(v.w) }; }
inline constexpr vec4 min(const vec4& lhs, const vec4& rhs) { return { min(lhs.x, rhs.x), min(lhs.y, rhs.y), min(lhs.z, rhs.z), min(lhs.w, rhs.w) }; }
inline constexpr vec4 max(const vec4& lhs, const vec4& rhs) { return { max(lhs.x, rhs.x), max(lhs.y, rhs.y), max(lhs.z, rhs.z), max(lhs.w, rhs.w) }; }
inline constexpr float min(const vec4& v) { return min(v.x, min(v.y, min(v.z, v.w))); }
inline constexpr float max(const vec4& v) { return max(v.x, max(v.y, max(v.z, v.w))); }
inline constexpr vec4 clamp(const vec4& v, const vec4& lo, const vec4& hi) { return { clamp(v.x, lo.x, hi.x), clamp(v.y, lo.y, hi.y), clamp(v.z, lo.z, hi.z), clamp(v.w, lo.w, hi.w) }; }
inline           float uangle(const vec4& a, const vec4& b) { float d = dot(a, b); return d > 1 ? 0 : std::acos(d < -1 ? -1 : d); }
inline constexpr vec4 lerp(const vec4& a, const vec4& b, float t) { return a + (b - a) * t; }
inline constexpr vec4 lerp(const vec4& a, const vec4& b, const vec4& t) { return { lerp(a.x, b.x, t.x), lerp(a.y, b.y, t.y), lerp(a.z, b.z, t.z), lerp(a.w, b.w, t.w) }; }
inline           vec4 nlerp(const vec4& a, const vec4& b, float t) { return normalize(lerp(a, b, t)); }
inline           vec4 slerp(const vec4& a, const vec4& b, float t) { float th = uangle(a, b); return th == 0 ? a : a * (std::sin(th * (1 - t)) / std::sin(th)) + b * (std::sin(th * t) / std::sin(th)); }


typedef vec4 Color;
typedef vec4 quat;

inline constexpr quat qconj(const quat& q) { return { -q.x, -q.y, -q.z, q.w }; }
inline constexpr quat qinv(const quat& q) { return qconj(q) / length2(q); }
inline constexpr quat qmul(const quat& a, const quat& b) { return {a.x*b.w+a.w*b.x+a.y*b.z-a.z*b.y, a.y*b.w+a.w*b.y+a.z*b.x-a.x*b.z, a.z*b.w+a.w*b.z+a.x*b.y-a.y*b.x, a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z}; }
inline constexpr vec3 qxdir (const quat& q) { return {q.w*q.w+q.x*q.x-q.y*q.y-q.z*q.z, (q.x*q.y+q.z*q.w)*2, (q.z*q.x-q.y*q.w)*2}; }
inline constexpr vec3 qydir (const quat& q) { return {(q.x*q.y-q.z*q.w)*2, q.w*q.w-q.x*q.x+q.y*q.y-q.z*q.z, (q.y*q.z+q.x*q.w)*2}; }
inline constexpr vec3 qzdir (const quat& q) { return {(q.z*q.x+q.y*q.w)*2, (q.y*q.z-q.x*q.w)*2, q.w*q.w-q.x*q.x-q.y*q.y+q.z*q.z}; }
inline constexpr vec3 qrot (const quat& q, const vec3& v) { return qxdir(q)*v.x + qydir(q)*v.y + qzdir(q)*v.z; }
inline           quat qnlerp(const quat& a, const quat& b, float t) { return nlerp(a, dot(a, b) < 0 ? -b : b, t); }
inline           quat qslerp(const quat& a, const quat& b, float t) { return slerp(a, dot(a, b) < 0 ? -b : b, t); }
inline           float qangle(const quat& q) { return std::atan2(length(q.xyz()), q.w)*2; }
inline           vec3 qaxis (const quat& q) { return normalize(q.xyz()); }
inline           quat qaxisangle(const vec3 & axis, float radianAngle) {
	// Assumes normalized axis
	const float halfAngle = radianAngle * 0.5f, sin_a = std::sin(halfAngle);
	quat q(axis.x * sin_a, axis.y * sin_a, axis.z * sin_a, std::cos(halfAngle));
	return normalize(q);
}
inline quat qeuler(const vec3& eulerDegrees) {
	const vec3 halfAngles = eulerDegrees * DEG_TO_RAD * 0.5f;
	float cx = cos(halfAngles.x);
	float sx = sin(halfAngles.x);
	float cy = cos(halfAngles.y);
	float sy = sin(halfAngles.y);
	float cz = cos(halfAngles.z);
	float sz = sin(halfAngles.z);
	return {
		sx * cy * cz - cx * sy * sz,
		cx * sy * cz + sx * cy * sz,
		cx * cy * sz - sx * sy * cz,
		cx * cy * cz + sx * sy * sz
	};
}
inline vec3 qeuler(const quat& q) { // returns euler degrees
	vec3 ret;
	// roll
	float sinr_cosp = 2.f * (q.w * q.x + q.y * q.z);
	float cosr_cosp = 1.f - 2.f * (q.x * q.x + q.y * q.y);
	ret.x = std::atan2(sinr_cosp, cosr_cosp);
	// pitch
	float sinp = 2.f * (q.w * q.y - q.z * q.x);
	if (std::abs(sinp) >= 1.f)
		ret.y = std::copysign(PI / 2.f, sinp); // use 90 degrees if out of range
	else
		ret.y = std::asin(sinp);
	// yaw
	float siny_cosp = 2.f * (q.w * q.z + q.x * q.y);
	float cosy_cosp = 1.f - 2.f * (q.y * q.y + q.z * q.z);
	ret.z = std::atan2(siny_cosp, cosy_cosp);
	return ret * RAD_TO_DEG;
}
inline quat qfromto(vec3 from, vec3 to) { // assumes normalized
	float cosTheta = dot(from, to) + 1.f;
	if (cosTheta < 0.0001f) {
		return abs(from.x) > abs(from.z)
			? normalize({ -from.x, from.y, 0, 0 })
			: normalize({ 0, -from.z, from.y, 0 });
	}
	vec3 axis = cross(from, to);
	return normalize({ axis.x, axis.y, axis.z, cosTheta });
}
inline quat qlookat(const vec3& dir, const vec3& desiredUp = { 0, 1, 0 }, const vec3& objectForward = { 0, 0, 1 }) { // assumes normalized
	vec3 up = desiredUp;
	vec3 right = orthonormalize(dir, up);

	quat rot1 = qfromto(objectForward, dir);
	vec3 newUp = qrot(rot1, desiredUp);
	quat rot2 = qfromto(newUp, up);
	return qmul(rot2, rot1);
}

struct Bounds {
	vec3 min;
	vec3 max;

	constexpr vec3 center() const { return (min + max) * 0.5f; }
	constexpr vec3 size() const { return max - min; }
	constexpr vec3 extents() const { return (max - min) * 0.5f; }
	constexpr Bounds& expand(const vec3& amount ) { min -= amount; max += amount; return *this; }
	constexpr bool contains(const vec3& p) const {
		return p.x >= min.x && p.y >= min.y && p.z >= min.z && p.x <= max.x && p.y <= max.y && p.z <= max.z;
	}
};


struct TriFace
{
	union {
		struct { uint a, b, c; };
		uint v[3];
	};
	constexpr TriFace(uint a_, uint b_, uint c_): a(a_), b(b_), c(c_) {}
	constexpr void offset(int offset) { a += offset; b += offset; c += offset; }
};

struct QuadFace
{
	union {
		struct { uint a, b, c, d; };
		uint v[4];
	};
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
	std::vector<Color> colors;

	uint numIndices() const { return triangles.empty() ? quads.size() * 4 : triangles.size() * 3; }
	uint* indices() { return triangles.empty() ? (uint*)quads.data() : (uint*)triangles.data(); }
	const uint* indices() const { return triangles.empty() ? (uint*)quads.data() : (uint*)triangles.data(); }
};

void bounds(const Geometry& geo, vec3& lower, vec3& upper);
Bounds bounds(const Geometry& geo);
float boundingSphere(const Geometry& geo);

Geometry& scale(Geometry& geo, vec3 s);
Geometry& translate(Geometry& geo, vec3 v);
Geometry& center(Geometry& geo);
Geometry& rotate(Geometry& geo, vec3 eulerDegrees);
Geometry& rotate(Geometry& geo, quat rot);
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
inline Geometry& rotatex(Geometry& geo, float degrees) { return rotate(geo, {degrees, 0, 0}); }
inline Geometry& rotatey(Geometry& geo, float degrees) { return rotate(geo, {0, degrees, 0}); }
inline Geometry& rotatez(Geometry& geo, float degrees) { return rotate(geo, {0, 0, degrees}); }

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
Geometry plane(vec2 size, vec3 normal);
// Cylinder
Geometry cylinder(float height, float radius, vec3 direction = {0, 1, 0}, int segments = 16);
Geometry cylinder(vec3 start, vec3 end, float radius, int segments = 16);
// Cone
Geometry cone(float height, float baseRadius, float tipRadius = 0.f, vec3 direction = {0, 1, 0}, int segments = 16);
Geometry cone(vec3 base, vec3 tip, float baseRadius, float tipRadius = 0.f, int segments = 16);

Geometry heightmap(float* heights, int w, int h, vec3 scale = vec3(1), bool center = true);
Geometry uvMesh(const Geometry& geo);
Geometry spin(const vec3* path, int numPathPoints, int spinSteps, const vec3& axis = {0, 1, 0}, bool autoMerge = false);
inline Geometry spin(const std::vector<vec3>& path, int spinSteps, const vec3& axis = {0, 1, 0}, bool autoMerge = false) {
	return spin(path.data(), path.size(), spinSteps, axis, autoMerge);
}

Geometry& subdivide(Geometry& geo, int amount = 1);
Geometry& weld(Geometry& geo, float maxNormalAngleDeg = 30.f, double tolerance = 0.0001);
Geometry& unweld(Geometry& geo, bool faceNormals);
Geometry& removeDuplicateFaces(Geometry& geo); // Only concerns positions, not normals/winding/texcoords etc.
Geometry& removeFacesWithNormal(Geometry& geo, vec3 normal, float maxAngleDegrees);
Geometry& removeFacesBehindPlane(Geometry& geo, vec3 planePoint, vec3 planeNormal);
Geometry& removeOrphanVertices(Geometry& geo);
Geometry& triangulate(Geometry& geo);
Geometry& displaceAlongNormals(Geometry& geo, float minAmount, float maxAmount, int seed = 0);
Geometry& displacementMap(Geometry& geo, float* displacementTexture, int w, int h, float height = 1.0f);
Geometry& displacementMap(Geometry& geo, std::function<float(vec2)> sampler, float height = 1.0f);
Geometry& bakeTextureToVertexColors(Geometry& geo, std::function<Color(vec2)> sampler);
Geometry& assignConstantVertexColor(Geometry& geo, Color color);
Geometry& calculateNormals(Geometry& geo);
Geometry& normalizeNormals(Geometry& geo);
Geometry& flipNormals(Geometry& geo);
Geometry& generateBoxTexcoords(Geometry& geo);

enum class MeshFormatStyle { ASCII, Binary };

bool writeObj(const Geometry& geo, const char* filename = "out.obj", bool optimize = true);
bool writeStl(const Geometry& geo, const char* filename = "out.stl", MeshFormatStyle formatStyle = MeshFormatStyle::Binary);

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
		p.x -= clamp(p.x, 0.f, h);
		return length(p) - r;
	}
	inline float capsuley(vec3 p, float h, float r) {
		p.y -= clamp(p.y, 0.f, h);
		return length(p) - r;
	}
	inline float capsulez(vec3 p, float h, float r) {
		p.z -= clamp(p.z, 0.f, h);
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
		vec2 d = abs(vec2(length(vec2(p.y, p.z)), p.x)) - vec2(r, h);
		return min(max(d.x, d.y), 0.f) + length(max(d, 0.f));
	}
	inline float cylindery(vec3 p, float h, float r) {
		vec2 d = abs(vec2(length(vec2(p.x, p.z)), p.y)) - vec2(r, h);
		return min(max(d.x, d.y), 0.f) + length(max(d, 0.f));
	}
	inline float cylinderz(vec3 p, float h, float r) {
		vec2 d = abs(vec2(length(vec2(p.x, p.y)), p.z)) - vec2(r, h);
		return min(max(d.x, d.y), 0.f) + length(max(d, 0.f));
	}

	inline constexpr vec3 translate(vec3 p, vec3 amount) { return p - amount; }
	inline constexpr float opOnion(float sdf, float thickness) { return abs(sdf) - thickness; }
	inline constexpr float opRound(float sdf, float r) { return sdf - r; }
	inline constexpr float opUnion(float a, float b) { return min(a, b); }
	inline constexpr float opSubtract(float a, float b) { return max(a, -b); }
	inline constexpr float opIntersect(float a, float b) { return max(a, b); }
	inline constexpr float opSmoothUnion(float a, float b, float k) {
		const float h = saturate(0.5f + 0.5f * (b - a) / k);
		return lerp(b, a, h) - k * h * (1.f - h);
	}
	inline constexpr float opSmoothSubtraction(float a, float b, float k) {
		const float h = saturate(0.5f - 0.5f * (b + a) / k);
		return lerp(-b, a, h) + k * h * (1.f - h);
	}
	inline constexpr float opSmoothIntersection(float a, float b, float k) {
		const float h = saturate(0.5f - 0.5f * (b - a) / k);
		return lerp(b, a, h) + k * h * (1.f - h);
	}

	struct SDF {
		enum SDFInstruction : unsigned char {
			NOOP,
			PRIM_SPHERE,
			PRIM_BOX,
			PRIM_PLANE,
			PRIM_CAPSULE,
			PRIM_CAPSULEX,
			PRIM_CAPSULEY,
			PRIM_CAPSULEZ,
			PRIM_CYLINDER,
			PRIM_CYLINDERX,
			PRIM_CYLINDERY,
			PRIM_CYLINDERZ,

			OP_UNION,
			OP_SUBTRACT,
			OP_INTERSECT,
			OP_UNION_SMOOTH,
			OP_SUBTRACT_SMOOTH,
			OP_INTERSECT_SMOOTH,
			OP_ROUND,
			OP_ONION,
			OP_SYMMETRY,

			PUSH_POS,
			POP_POS
		};

		SDF& combine(const SDF& other, SDFInstruction op = SDF::OP_UNION);

		float operator()(vec3 pos);
		bool empty() const { return code.empty(); }

		SDF& round(float r);
		SDF& onion(float thinckess);
		SDF& symmetry(bool x, bool y, bool z);
		SDF& sphere(vec3 pos, float r);
		SDF& box(vec3 pos, vec3 bounds);
		SDF& plane(vec3 pos, vec3 n, float d);
		SDF& capsule(vec3 a, vec3 b, float r);
		SDF& capsulex(vec3 pos, float h, float r);
		SDF& capsuley(vec3 pos, float h, float r);
		SDF& capsulez(vec3 pos, float h, float r);
		SDF& cylinder(vec3 a, vec3 b, float r);
		SDF& cylinderx(vec3 pos, float h, float r);
		SDF& cylindery(vec3 pos, float h, float r);
		SDF& cylinderz(vec3 pos, float h, float r);

	private:
		// Note: all these need to be handled in combine()
		std::vector<SDFInstruction> code;
		std::vector<float> params;
		std::vector<Geometry> geos;
		std::vector<Bounds> geoBounds;
	};

	typedef std::function<float(vec3)> SDFFunction;
	Geometry polygonize(vec3 start, vec3 end, vec3 step, SDFFunction sdf);
	std::vector<float> slice(vec2 size, vec3 topLeft, vec3 topRight, vec3 bottomLeft, vec3 bottomRight, float offset, SDFFunction sdf);

} // namespace sdf


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
Geometry polygonize(const VoxelGrid& voxelGrid, float voxelSize = 1.0f);

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
	template<> struct hash<gengeo::vec4> {
		size_t operator()(const gengeo::vec4& v) const {
			size_t seed = 62339909;
			std::hash<float> hasher;
			gengeo::hash_combine(seed, hasher(v.x));
			gengeo::hash_combine(seed, hasher(v.y));
			gengeo::hash_combine(seed, hasher(v.z));
			gengeo::hash_combine(seed, hasher(v.w));
			return seed;
		}
	};
}