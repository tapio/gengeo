#include "gengeo.hpp"
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <ctime>
#include <cassert>
#include <random>
#include <numeric>
#include <unordered_map>
#include "../third-party/tinyformat.h"

namespace gengeo {

#include "marchingcubes.h"

enum class UvConvention { DX, OGL };
static constexpr UvConvention UvStyle = UvConvention::OGL;
static constexpr float FlipYForOGL(float y, UvConvention style) { return style == UvConvention::OGL ? 1.0f - y : y; }
static constexpr float FlipYForDX(float y, UvConvention style) { return style == UvConvention::DX ? 1.0f - y : y; }
static inline void addQuadUVs(std::vector<vec2>& texcoords, UvConvention style)
{
	const float top = style == UvConvention::DX ? 0.0f : 1.0f;
	const float bottom = style == UvConvention::DX ? 1.0f : 0.0f;
	texcoords.emplace_back(0.f, bottom);
	texcoords.emplace_back(1.f, bottom);
	texcoords.emplace_back(1.f, top);
	texcoords.emplace_back(0.f, top);
}
Geometry cube()
{
	static constexpr vec3 verts[8] = {
		{-0.5f, -0.5f, -0.5f}, // 0
		{-0.5f,  0.5f, -0.5f}, // 1
		{ 0.5f,  0.5f, -0.5f}, // 2
		{ 0.5f, -0.5f, -0.5f}, // 3
		{-0.5f, -0.5f,  0.5f}, // 4
		{-0.5f,  0.5f,  0.5f}, // 5
		{ 0.5f,  0.5f,  0.5f}, // 6
		{ 0.5f, -0.5f,  0.5f}  // 7
	};
	static constexpr QuadFace quads[6] = {
		{4, 7, 6, 5}, // front
		{3, 0, 1, 2}, // back
		{7, 3, 2, 6}, // right
		{5, 6, 2, 1}, // top
		{0, 4, 5, 1}, // left
		{7, 4, 0, 3}  // bottom
	};
	static constexpr vec3 normals[6] = {
		{ 0,  0,  1}, // front
		{ 0,  0, -1}, // back
		{ 1,  0,  0}, // right
		{ 0,  1,  0}, // top
		{-1,  0,  0}, // left
		{ 0, -1,  0}  // bottom
	};
	Geometry geo;
	geo.quads.reserve(6);
	geo.positions.reserve(4*6);
	geo.texcoords.reserve(4*6);
	geo.normals.reserve(4*6);
	for (int i = 0; i < 6; i++) {
		geo.quads.emplace_back(i*4+0, i*4+1, i*4+2, i*4+3);
		geo.positions.emplace_back(verts[quads[i].a]);
		geo.positions.emplace_back(verts[quads[i].b]);
		geo.positions.emplace_back(verts[quads[i].c]);
		geo.positions.emplace_back(verts[quads[i].d]);
		geo.normals.emplace_back(normals[i]);
		geo.normals.emplace_back(normals[i]);
		geo.normals.emplace_back(normals[i]);
		geo.normals.emplace_back(normals[i]);
		addQuadUVs(geo.texcoords, UvStyle);
	}
	return geo;
}

Geometry sphere(int subdivs)
{
	Geometry geo = cube();
	subdivide(scale(geo, 2.f), subdivs);
	for (uint i = 0; i < geo.positions.size(); ++i) {
		vec3& pos = geo.positions[i];
		float x = pos.x, y = pos.y, z = pos.z;
		// http://mathproofs.blogspot.co.uk/2005/07/mapping-cube-to-sphere.html
		pos.x = x * sqrt(1.f - y*y*0.5f - z*z*0.5f + y*y*z*z/3.f);
		pos.y = y * sqrt(1.f - z*z*0.5f - x*x*0.5f + z*z*x*x/3.f);
		pos.z = z * sqrt(1.f - x*x*0.5f - y*y*0.5f + x*x*y*y/3.f);
		// geo.positions[i] = normalize(pos); // simple normalization
		geo.normals[i] = normalize(pos);
	}
	return geo;
}

Geometry uvsphere(int latSlices, int longSlices)
{
	if (latSlices < 3)
		latSlices = 3;
	if (longSlices < 3)
		longSlices = 3;

	std::vector<vec3> profile; // TODO: Would be nice to do on stack
	float latAngleStep = PI / (latSlices - 1);
	for (int i = 0; i < latSlices; i++) {
		float x = std::sin(i * latAngleStep);
		float y = -std::cos(i * latAngleStep);
		profile.emplace_back(x, y, 0.f);
	}

	constexpr vec3 up = {0, 1, 0}; // TODO: configure
	return spin(profile.data(), profile.size(), longSlices, up, false);
}

Geometry circley(int points)
{
	if (points < 3)
		points = 3;
	int total = points + 1;
	Geometry geo;
	geo.positions.reserve(total);
	geo.texcoords.reserve(total);
	geo.normals.reserve(total);
	geo.triangles.reserve(points);
	geo.positions.emplace_back(0.f, 0.f, 0.f);
	geo.texcoords.emplace_back(0.5f, 0.5f);
	float angleStep = TWOPI / points;
	for (int i = 0; i < points; i++) {
		geo.triangles.emplace_back(0, ((i+1) % points) + 1, i+1);
		float x = std::cos(i * angleStep);
		float y = std::sin(i * angleStep);
		float u = (x + 1) * 0.5f;
		float v = 1.0f - ((y + 1) * 0.5f);
		geo.positions.emplace_back(x, 0.f, y);
		geo.texcoords.emplace_back(u, v);
	}
	geo.normals.assign(total, {0, 1, 0});
	return geo;
}

Geometry planexy()
{
	Geometry geo;
	geo.quads.emplace_back(0, 1, 2, 3);
	geo.positions.emplace_back(-0.5f, -0.5f, 0.f);
	geo.positions.emplace_back( 0.5f, -0.5f, 0.f);
	geo.positions.emplace_back( 0.5f,  0.5f, 0.f);
	geo.positions.emplace_back(-0.5f,  0.5f, 0.f);
	addQuadUVs(geo.texcoords, UvStyle);
	geo.normals.assign(4, { 0.f, 0.f, 1.f });
	return geo;
}

Geometry planexz()
{
	Geometry geo;
	geo.quads.emplace_back(0, 1, 2, 3);
	geo.positions.emplace_back(-0.5f, 0.f,  0.5f);
	geo.positions.emplace_back( 0.5f, 0.f,  0.5f);
	geo.positions.emplace_back( 0.5f, 0.f, -0.5f);
	geo.positions.emplace_back(-0.5f, 0.f, -0.5f);
	addQuadUVs(geo.texcoords, UvStyle);
	geo.normals.assign(4, { 0.f, 1.f, 0.f });
	return geo;
}

Geometry planeyz()
{
	Geometry geo;
	geo.quads.emplace_back(0, 1, 2, 3);
	geo.positions.emplace_back(0.f, -0.5f,  0.5f);
	geo.positions.emplace_back(0.f, -0.5f, -0.5f);
	geo.positions.emplace_back(0.f,  0.5f, -0.5f);
	geo.positions.emplace_back(0.f,  0.5f,  0.5f);
	addQuadUVs(geo.texcoords, UvStyle);
	geo.normals.assign(4, { 1.f, 0.f, 0.f });
	return geo;
}

Geometry plane(vec2 size, vec3 normal)
{
	if (normal == vec3(0, 0, 1)) {
		Geometry geo = planexy();
		return scale(geo, vec3(size.x, size.y, 1.f));
	} else if (normal == vec3(0, 1, 0)) {
		Geometry geo = planexz();
		return scale(geo, vec3(size.x, 1.f, size.y));
	} else if (normal == vec3(1, 0, 0)) {
		Geometry geo = planeyz();
		return scale(geo, vec3(1.f, size.x, size.y));
	}
	Geometry geo = planexy();
	geo = scale(geo, vec3(size.x, size.y, 1.f));
	//quat rot = qfromto(vec3(0, 0, 1), normalize(normal));
	quat rot = qlookat(normalize(normal));
	return rotate(geo, rot);
}

Geometry cylinder(float height, float radius, vec3 direction, int segments)
{
	direction = normalize(direction);
	const vec3 hh = direction * (height * 0.5f);
	vec3 someAxis(0, 0, 1);
	if (abs(dot(direction, someAxis)) > 0.95f) // Almost parallel, let's pick another axis
		someAxis = {0, 1, 0};
	const vec3 side = normalize(cross(direction, someAxis)) * radius;
	const vec3 profile[] = { -hh, -hh + side, hh + side, hh };
	return spin(profile, 4, segments, direction, false);
}

Geometry cylinder(vec3 start, vec3 end, float radius, int segments)
{
	vec3 direction = normalize(end - start);
	vec3 someAxis(0, 0, 1);
	if (abs(dot(direction, someAxis)) > 0.95f) // Almost parallel, let's pick another axis
		someAxis = {0, 1, 0};
	const vec3 side = normalize(cross(direction, someAxis)) * radius;
	const vec3 profile[] = { start, start + side, end + side, end };
	return spin(profile, 4, segments, direction, false);
}

Geometry cone(float height, float baseRadius, float tipRadius, vec3 direction, int segments)
{
	return cone({}, normalize(direction) * height, baseRadius, tipRadius, segments);
}

Geometry cone(vec3 base, vec3 tip, float baseRadius, float tipRadius, int segments)
{
	vec3 direction = normalize(tip - base);
	vec3 someAxis(0, 0, 1);
	if (abs(dot(direction, someAxis)) > 0.95f) // Almost parallel, let's pick another axis
		someAxis = {0, 1, 0};
	const vec3 side = normalize(cross(direction, someAxis));
	if (tipRadius == 0.f) {
		const vec3 profile[] = { base, base + side * baseRadius, tip };
		return spin(profile, 3, segments, direction, false);
	} else {
		const vec3 profile[] = { base, base + side * baseRadius, tip + side * tipRadius, tip };
		return spin(profile, 4, segments, direction, false);
	}
}

Geometry heightmap(float* heights, int w, int h, vec3 scale, bool center)
{
	int numVerts = w * h;
	int lasti = w - 1;
	int lastj = h - 1;
	float offsetx = center ? w * -0.5f : 0.f;
	float offsetz = center ? h * -0.5f : 0.f;
	Geometry geo;
	geo.positions.resize(numVerts);
	geo.texcoords.resize(numVerts);
	geo.normals.resize(numVerts);
	geo.quads.reserve(numVerts);
	for (int j = 0; j < h; ++j) {
		for (int i = 0; i < w; ++i) {
			// Create the vertex
			float x = i + offsetx;
			float z = j + offsetz;
			float y = heights[j * w + i];
			int vert = j * w + i;
			geo.positions[vert] = vec3(x, y, z) * scale;
			geo.texcoords[vert] = vec2((float)i / lasti, (float)j / lastj);
			geo.normals[vert] = vec3(0, 1, 0);
			// Indexed faces
			if (i == lasti || j == lastj)
				continue;
			uint a = i + w * j;
			uint b = i + w * (j + 1);
			uint c = (i + 1) + w * (j + 1);
			uint d = (i + 1) + w * j;
			//uint triangles[] = { a, b, d, b, c, d };
			//geo.triangles.insert(geo.triangles.end(), triangles, triangles + 6);
			geo.quads.push_back({ a, b, c, d });
		}
	}
	calculateNormals(geo);
	return geo;
}

Geometry uvMesh(const Geometry& geo)
{
	if (geo.texcoords.empty())
		return Geometry();

	Geometry ret;
	bool quads = geo.quads.size() > 0;
	uint numIndices = quads ? geo.quads.size() * 4 : geo.triangles.size() * 3;
	uint* indices = quads ? (uint*)geo.quads.data() : (uint*)geo.triangles.data();
	ret.positions.reserve(numIndices);
	ret.texcoords.reserve(numIndices);
	ret.normals.reserve(numIndices);
	// TODO: Fix winding for half of them...
	for (uint i = 0; i < numIndices; ++i) {
		vec2 uv = geo.texcoords[indices[i]];
		ret.positions.emplace_back(uv.x, uv.y, 0.f);
		ret.texcoords.emplace_back(uv.x, uv.y);
		ret.normals.emplace_back(0, 0, 1);
	}
	if (quads) {
		ret.quads.resize(geo.quads.size(), {0, 0, 0, 0});
		std::iota((uint*)ret.quads.data(), (uint*)ret.quads.data() + numIndices, 0); // Fill indices with ascending number
	} else {
		ret.triangles.resize(geo.triangles.size(), {0, 0, 0});
		std::iota((uint*)ret.triangles.data(), (uint*)ret.triangles.data() + numIndices, 0); // Fill indices with ascending number
	}
	return ret;
}

Geometry spin(const vec3* path, int numPathPoints, int steps, const vec3& axis, bool autoMerge)
{
	Geometry ret;
	if (numPathPoints < 2 || steps < 2)
		return ret;
	float angle = TWOPI / steps;

	// TODO: Detect up/down path direction for UVs and faces...
	// TODO: Triangles would be better shapes that touch the axis
	// Currently assumes from bottom to up path

	ret.positions.assign(path, path + numPathPoints);
	ret.positions.reserve(numPathPoints * steps);
	ret.quads.reserve((numPathPoints - 1) * steps);
	ret.texcoords.reserve(ret.positions.size());

	// Texcoord precompute
	ret.texcoords.emplace_back(0, 0);
	float pathLength = 0.f;
	for (int i = 1; i < numPathPoints; ++i) {
		const vec3 prev = path[i - 1];
		const vec3 cur = path[i];
		pathLength += distance(prev, cur);
		ret.texcoords.emplace_back(0, pathLength);
	}
	pathLength = 1.f / pathLength;
	for (vec2& uv : ret.texcoords)
		uv.y = FlipYForDX(uv.y * pathLength, UvStyle);

	// Vertices
	for (int i = 1; i < steps; ++i) {
		quat q = qaxisangle(axis, angle * i);
		for (int j = 0; j < numPathPoints; ++j) {
			ret.positions.push_back(qrot(q, path[j]));
			ret.texcoords.emplace_back((float)i / steps, ret.texcoords[j].y);
		}
	}
	if (!autoMerge) { // Last slice duplicated if not merging
		for (int i = 0; i < numPathPoints; ++i) {
			ret.positions.push_back(path[i]);
			ret.texcoords.emplace_back(1, ret.texcoords[i].y);
		}
	}

	// Faces (quads)
	int modOp = autoMerge ? steps : (steps + 1);
	for (int step = 1; step <= steps; ++step) {
		int curOffset = (step % modOp) * numPathPoints;
		int prevOffset = (step - 1) * numPathPoints;
		for (int index = 1; index < numPathPoints; ++index) {
			//ret.quads.emplace_back(prevOffset + index - 1, prevOffset + index, curOffset + index, curOffset + index - 1);
			ret.quads.emplace_back(prevOffset + index - 1, curOffset + index - 1, curOffset + index, prevOffset + index);
		}
	}

	return calculateNormals(ret);
}



void bounds(const Geometry& geo, vec3& lower, vec3& upper)
{
	vec3 lbound(std::numeric_limits<float>::max());
	vec3 ubound(-std::numeric_limits<float>::max());
	for (auto pos : geo.positions) {
		lbound = min(pos, lbound);
		ubound = max(pos, ubound);
	}
	lower = lbound;
	upper = ubound;
}

Bounds bounds(const Geometry& geo)
{
	Bounds b;
	bounds(geo, b.min, b.max);
	return b;
}

float boundingSphere(const Geometry& geo)
{
	float maxRadiusSq = 0;
	for (auto pos : geo.positions)
		maxRadiusSq = max(maxRadiusSq, length2(pos));
	return std::sqrt(maxRadiusSq);
}

Geometry& scale(Geometry& geo, vec3 s)
{
	if (s == vec3(1.f))
		return geo;
	// TODO: Normals (non-uniform, mirroring)
	for (auto& pos : geo.positions)
		pos *= s;
	return geo;
}

Geometry& translate(Geometry& geo, vec3 v)
{
	if (v == vec3(0.f))
		return geo;
	for (auto& pos : geo.positions)
		pos += v;
	return geo;
}

Geometry& center(Geometry& geo)
{
	vec3 lower, upper;
	bounds(geo, lower, upper);
	translate(geo, -(lower + upper) * 0.5f);
	return geo;
}

Geometry& rotate(Geometry& geo, vec3 eulerDegrees)
{
	if (eulerDegrees == vec3(0.f))
		return geo;
	return rotate(geo, qeuler(eulerDegrees));
}

Geometry& rotate(Geometry& geo, quat rot)
{
	if (rot == quat(0.f, 0.f, 0.f, 1.f))
		return geo;
	for (auto& pos : geo.positions)
		pos = qrot(rot, pos);
	for (auto& n : geo.normals)
		n = qrot(rot, n);
	return geo;
}

Geometry& concat(Geometry& geo, const Geometry& other)
{
	uint indexOffset = geo.positions.size();
	geo.positions.insert(geo.positions.end(), other.positions.begin(), other.positions.end());
	geo.texcoords.insert(geo.texcoords.end(), other.texcoords.begin(), other.texcoords.end());
	geo.normals.insert(geo.normals.end(), other.normals.begin(), other.normals.end());
	geo.colors.insert(geo.colors.end(), other.colors.begin(), other.colors.end());

	uint firstNewQuad = geo.quads.size();
	geo.quads.insert(geo.quads.end(), other.quads.begin(), other.quads.end());
	for (uint i = firstNewQuad; i < geo.quads.size(); ++i)
		geo.quads[i].offset(indexOffset);

	uint firstNewTri = geo.triangles.size();
	geo.triangles.insert(geo.triangles.end(), other.triangles.begin(), other.triangles.end());
	for (uint i = firstNewTri; i < geo.triangles.size(); ++i)
		geo.triangles[i].offset(indexOffset);
	return geo;
}

template<typename T>
static void subdiv(const QuadFace quad, T& destination, const T& source)
{
	if (source.size() <= quad.a)
		return;
	/*0*/ destination.emplace_back(source[quad.a]);
	/*1*/ destination.emplace_back((source[quad.a] + source[quad.b]) * 0.5f);
	/*2*/ destination.emplace_back(source[quad.b]);
	/*3*/ destination.emplace_back((source[quad.b] + source[quad.c]) * 0.5f);
	/*4*/ destination.emplace_back(source[quad.c]);
	/*5*/ destination.emplace_back((source[quad.c] + source[quad.d]) * 0.5f);
	/*6*/ destination.emplace_back(source[quad.d]);
	/*7*/ destination.emplace_back((source[quad.d] + source[quad.a]) * 0.5f);
	/*8*/ destination.emplace_back((source[quad.a] + source[quad.b] + source[quad.c] + source[quad.d]) * 0.25f);
}

Geometry& subdivide(Geometry& geo, int amount)
{
	while (amount--) {
		Geometry subdivided;
		subdivided.positions.reserve(geo.positions.size() * 4);
		subdivided.texcoords.reserve(geo.texcoords.size() * 4);
		subdivided.normals.reserve(geo.normals.size() * 4);
		subdivided.colors.reserve(geo.colors.size() * 4);
		subdivided.quads.reserve(geo.quads.size() * 4);
		for (const auto& quad : geo.quads) {
			uint base = subdivided.positions.size();
			subdiv(quad, subdivided.positions, geo.positions);
			subdiv(quad, subdivided.texcoords, geo.texcoords);
			subdiv(quad, subdivided.normals, geo.normals);
			subdiv(quad, subdivided.colors, geo.colors);
			subdivided.quads.emplace_back(base + 0, base + 1, base + 8, base + 7);
			subdivided.quads.emplace_back(base + 1, base + 2, base + 3, base + 8);
			subdivided.quads.emplace_back(base + 8, base + 3, base + 4, base + 5);
			subdivided.quads.emplace_back(base + 7, base + 8, base + 5, base + 6);
		}
		geo = std::move(subdivided);
		normalizeNormals(geo);
	}
	return weld(geo);
}

static void copyVertex(Geometry& dst, Geometry& src, uint srcIndex)
{
	dst.positions.emplace_back(src.positions[srcIndex]);
	if (srcIndex < src.texcoords.size())
		dst.texcoords.emplace_back(src.texcoords[srcIndex]);
	if (srcIndex < src.normals.size())
		dst.normals.emplace_back(src.normals[srcIndex]);
	if (srcIndex < src.colors.size())
		dst.colors.emplace_back(src.colors[srcIndex]);
}

static vec3 roundToTolerance(vec3 pos, double tolerance)
{
	return {
		(float)(std::round((double)pos.x / tolerance) * tolerance),
		(float)(std::round((double)pos.y / tolerance) * tolerance),
		(float)(std::round((double)pos.z / tolerance) * tolerance)
	};
}

Geometry& weld(Geometry& geo, float maxNormalAngleDeg, double tolerance)
{
	// Round to tolerance
	// TODO: This can modify already welded vertices which is not desirable, use copies for hashing
	for (auto& pos : geo.positions)
		pos = roundToTolerance(pos, tolerance);
	float cosTheta = std::cos(maxNormalAngleDeg * DEG_TO_RAD);
	std::unordered_map<vec3, uint> hashmap;
	Geometry welded;
	welded.positions.reserve(geo.positions.size());
	welded.texcoords.reserve(geo.texcoords.size());
	welded.normals.reserve(geo.normals.size());
	welded.colors.reserve(geo.colors.size());

	auto handleVertex = [&, cosTheta](uint oldIndex) {
		uint toCreateIndex = welded.positions.size();
		uint newIndex = hashmap.try_emplace(geo.positions[oldIndex], toCreateIndex).first->second;
		if (newIndex < toCreateIndex && !geo.normals.empty() && dot(geo.normals[oldIndex], welded.normals[newIndex]) < cosTheta)
			newIndex = toCreateIndex;
		if (newIndex == toCreateIndex)
			copyVertex(welded, geo, oldIndex);
		return newIndex;
	};

	if (!geo.triangles.empty()) {
		for (TriFace oldFace : geo.triangles) {
			TriFace newFace(handleVertex(oldFace.a), handleVertex(oldFace.b), handleVertex(oldFace.c));
			welded.triangles.emplace_back(newFace);
		}
	} else {
		for (QuadFace oldFace : geo.quads) {
			QuadFace newFace(handleVertex(oldFace.a), handleVertex(oldFace.b), handleVertex(oldFace.c), handleVertex(oldFace.d));
			welded.quads.emplace_back(newFace);
		}
	}

	geo = std::move(welded);
	return geo;
}

Geometry& unweld(Geometry& geo, bool faceNormals)
{
	Geometry unwelded;
	unwelded.positions.reserve(geo.positions.size());
	unwelded.texcoords.reserve(geo.texcoords.size());
	unwelded.normals.reserve(geo.normals.size());
	unwelded.colors.reserve(geo.colors.size());

	uint numIndices = geo.numIndices();
	uint* indices = geo.indices();
	for (uint i = 0; i < numIndices; ++i) {
		uint index = indices[i];
		unwelded.positions.push_back(geo.positions[index]);
		if (index < geo.texcoords.size())
			unwelded.texcoords.push_back(geo.texcoords[index]);
		if (index < geo.normals.size())
			unwelded.normals.push_back(geo.normals[index]);
		if (index < geo.colors.size())
			unwelded.colors.push_back(geo.colors[index]);
	}

	if (!geo.triangles.empty()) {
		for (uint i = 0; i < geo.triangles.size(); ++i)
			unwelded.triangles.emplace_back(i*3, i*3 + 1, i*3 + 2);
	} else {
		for (uint i = 0; i < geo.quads.size(); ++i)
			unwelded.quads.emplace_back(i*4, i*4 + 1, i*4 + 2, i*4 + 3);
	}

	if (faceNormals)
		calculateNormals(unwelded);
	geo = std::move(unwelded);
	return geo;
}


static size_t hashTriFacePos(const Geometry& geo, const TriFace a) {
	constexpr double tolerance = 0.00001;
	size_t hasha = std::hash<vec3>{}(roundToTolerance(geo.positions[a.a], tolerance));
	size_t hashb = std::hash<vec3>{}(roundToTolerance(geo.positions[a.b], tolerance));
	size_t hashc = std::hash<vec3>{}(roundToTolerance(geo.positions[a.c], tolerance));
	// We sort the elements to make the hash not sensitive to vertex order
	if (hasha > hashc) std::swap(hasha, hashc);
	if (hasha > hashb) std::swap(hasha, hashb);
	if (hashb > hashc) std::swap(hashb, hashc);
	hash_combine(hasha, hashb);
	hash_combine(hasha, hashc);
	return hasha;
}

static size_t hashQuadFacePos(const Geometry& geo, const QuadFace a) {
	constexpr double tolerance = 0.00001;
	size_t hasha = std::hash<vec3>{}(roundToTolerance(geo.positions[a.a], tolerance));
	size_t hashb = std::hash<vec3>{}(roundToTolerance(geo.positions[a.b], tolerance));
	size_t hashc = std::hash<vec3>{}(roundToTolerance(geo.positions[a.c], tolerance));
	size_t hashd = std::hash<vec3>{}(roundToTolerance(geo.positions[a.d], tolerance));
	// We sort the elements to make the hash not sensitive to vertex order
	if (hasha > hashb) std::swap(hasha, hashb);
	if (hashc > hashd) std::swap(hashc, hashd);
	if (hasha > hashc) std::swap(hasha, hashc);
	if (hashb > hashd) std::swap(hashb, hashd);
	if (hashb > hashc) std::swap(hashb, hashc);
	hash_combine(hasha, hashb);
	hash_combine(hasha, hashc);
	hash_combine(hasha, hashd);
	return hasha;
}

Geometry& removeDuplicateFaces(Geometry& geo)
{
	if (!geo.triangles.empty()) {
		std::unordered_map<uint, bool> duplicates;
		std::unordered_map<uint, uint> faceHashes;
		for (uint i = 0; i < geo.triangles.size(); ++i) {
			size_t hash = hashTriFacePos(geo, geo.triangles[i]);
			auto insertion = faceHashes.try_emplace(hash, i);
			if (!insertion.second) {
				duplicates[i] = true;
				duplicates[insertion.first->second] = true;
			}
		}
		for (int i = geo.triangles.size() - 1; i >= 0; --i) {
			if (duplicates[i])
				geo.triangles.erase(geo.triangles.begin() + i);
		}
	} else if (!geo.quads.empty()) {
		std::unordered_map<uint, bool> duplicates;
		std::unordered_map<uint, uint> faceHashes;
		for (uint i = 0; i < geo.quads.size(); ++i) {
			size_t hash = hashQuadFacePos(geo, geo.quads[i]);
			auto insertion = faceHashes.try_emplace(hash, i);
			if (!insertion.second) {
				duplicates[i] = true;
				duplicates[insertion.first->second] = true;
			}
		}
		for (int i = geo.quads.size() - 1; i >= 0; --i) {
			if (duplicates[i])
				geo.quads.erase(geo.quads.begin() + i);
		}

	}
	return removeOrphanVertices(geo);
}

Geometry& removeFacesWithNormal(Geometry& geo, vec3 normal, float maxAngleDegrees)
{
	normal = normalize(normal);
	const float maxAngle = maxAngleDegrees * DEG_TO_RAD;
	if (!geo.triangles.empty()) {
		for (int i = geo.triangles.size() - 1; i >= 0; --i) {
			TriFace face = geo.triangles[i];
			vec3 n = triangleNormal(geo.positions[face.a], geo.positions[face.b], geo.positions[face.c]);
			if (uangle(n, normal) <= maxAngle)
				geo.triangles.erase(geo.triangles.begin() + i);
		}
	} else if (!geo.quads.empty()) {
		for (int i = geo.quads.size() - 1; i >= 0; --i) {
			QuadFace face = geo.quads[i];
			vec3 n = quadNormal(geo.positions[face.a], geo.positions[face.b], geo.positions[face.c], geo.positions[face.d]);
			if (uangle(n, normal) <= maxAngle)
				geo.quads.erase(geo.quads.begin() + i);
		}
	}
	return removeOrphanVertices(geo);
}

Geometry& removeFacesBehindPlane(Geometry& geo, vec3 planePoint, vec3 planeNormal)
{
	planeNormal = normalize(planeNormal);
	auto pointBehindPlane = [=](vec3 point) {
		return dot(planeNormal, normalize(planePoint - point)) > 0.f;
	};
	if (!geo.triangles.empty()) {
		for (int i = geo.triangles.size() - 1; i >= 0; --i) {
			TriFace face = geo.triangles[i];
			if (pointBehindPlane(geo.positions[face.a]) && pointBehindPlane(geo.positions[face.b]) && pointBehindPlane(geo.positions[face.c]))
				geo.triangles.erase(geo.triangles.begin() + i);
		}
	} else if (!geo.quads.empty()) {
		for (int i = geo.quads.size() - 1; i >= 0; --i) {
			QuadFace face = geo.quads[i];
			if (pointBehindPlane(geo.positions[face.a]) && pointBehindPlane(geo.positions[face.b]) && pointBehindPlane(geo.positions[face.c]) && pointBehindPlane(geo.positions[face.d]))
				geo.quads.erase(geo.quads.begin() + i);
		}
	}
	return removeOrphanVertices(geo);
}

Geometry& removeOrphanVertices(Geometry& geo)
{
	if (geo.positions.empty())
		return geo;
	std::unordered_map<uint, bool> inUse;
	uint numIndices = geo.numIndices();
	uint* indices = geo.indices();
	for (uint i = 0; i < numIndices; ++i)
		inUse[indices[i]] = true;
	for (int v = geo.positions.size() - 1; v >= 0; --v) {
		if (inUse[v])
			continue;
		// Remove vertex
		geo.positions.erase(geo.positions.begin() + v);
		if (v < geo.texcoords.size())
			geo.texcoords.erase(geo.texcoords.begin() + v);
		if (v < geo.normals.size())
			geo.normals.erase(geo.normals.begin() + v);
		if (v < geo.colors.size())
			geo.colors.erase(geo.colors.begin() + v);
		// Shift indices
		for (uint i = 0; i < numIndices; ++i)
			if (indices[i] >= (uint)v)
				--indices[i];
	}
	return geo;
}

Geometry& triangulate(Geometry& geo)
{
	if (!geo.quads.empty()) {
		geo.triangles.clear();
		geo.triangles.reserve(geo.quads.size() * 2);
		for (auto& face : geo.quads) {
			geo.triangles.emplace_back(face.a, face.b, face.c);
			geo.triangles.emplace_back(face.c, face.d, face.a);
		}
		geo.quads.clear();
	}
	return geo;
}

Geometry& displaceAlongNormals(Geometry& geo, float minAmount, float maxAmount, int seed)
{
	// Need to weld everything to avoid cracks...
	// TODO: Don't actually weld, just consider normals welded
	weld(geo, 180);
	if (geo.normals.empty())
		calculateNormals(geo);
	auto rand = std::bind(std::uniform_real_distribution<float>(0, 1), std::mt19937(seed ? seed : std::time(0)));
	auto& pos = geo.positions;
	auto& n = geo.normals;
	for (uint i = 0; i < pos.size(); i++)
		pos[i] += n[i] * lerp(minAmount, maxAmount, rand());
	return calculateNormals(geo);
}

static constexpr float sampleTexture(vec2 uv, float* tex, int w, int h)
{
	int x = int(uv.x * (w - 1)) % w;
	int y = int(uv.y * (h - 1)) % h;
	return tex[w * y + x];
}

Geometry& displacementMap(Geometry& geo, float* displacementTexture, int w, int h, float height)
{
	return displacementMap(geo, [displacementTexture, w, h](vec2 uv) {
		return sampleTexture(uv, displacementTexture, w, h);
	}, height);
}

Geometry& displacementMap(Geometry& geo, std::function<float(vec2)> sampler, float height)
{
	if (geo.texcoords.empty())
		return geo;
	// Need to weld everything to avoid cracks...
	// TODO: Don't actually weld, just consider normals welded
	weld(geo, 180);
	if (geo.normals.empty())
		calculateNormals(geo);
	auto& pos = geo.positions;
	auto& uv = geo.texcoords;
	auto& n = geo.normals;
	for (uint i = 0; i < pos.size(); i++)
		pos[i] += n[i] * sampler(uv[i]) * height;
	return calculateNormals(geo);
}

Geometry& bakeTextureToVertexColors(Geometry& geo, std::function<Color(vec2)> sampler)
{
	auto& pos = geo.positions;
	auto& color = geo.colors;
	auto& uv = geo.texcoords;
	if (uv.empty() || uv.size() != pos.size())
		return geo;
	color.resize(pos.size());
	for (uint i = 0; i < pos.size(); i++)
		color[i] = sampler(uv[i]);
	return geo;
}

Geometry& assignConstantVertexColor(Geometry& geo, Color color)
{
	geo.colors.assign(geo.positions.size(), color);
	return geo;
}

Geometry& calculateNormals(Geometry& geo)
{
	auto& triangles = geo.triangles;
	auto& quads = geo.quads;
	auto& positions = geo.positions;
	auto& normals = geo.normals;
	if (!triangles.empty()) {
		normals.clear();
		normals.resize(positions.size());
		for (auto tri : triangles) {
			vec3 normal = triangleNormal(positions[tri.a], positions[tri.b], positions[tri.c]);
			normals[tri.a] += normal;
			normals[tri.b] += normal;
			normals[tri.c] += normal;
		}
	} else if (!quads.empty()) {
		normals.clear();
		normals.resize(positions.size());
		for (auto quad : quads) {
			// TODO: Don't assume planar quads?
			#if 1
			vec3 normal = quadNormal(positions[quad.a], positions[quad.b], positions[quad.c], positions[quad.d]);
			normals[quad.a] += normal;
			normals[quad.b] += normal;
			normals[quad.c] += normal;
			normals[quad.d] += normal;
			#else
			normals[quad.a] += triangleNormal(positions[quad.a], positions[quad.b], positions[quad.d]);;
			normals[quad.b] += triangleNormal(positions[quad.b], positions[quad.a], positions[quad.c]);;
			normals[quad.c] += triangleNormal(positions[quad.c], positions[quad.b], positions[quad.d]);;
			normals[quad.d] += triangleNormal(positions[quad.d], positions[quad.a], positions[quad.c]);;
			#endif
		}
	}
	normalizeNormals(geo);
	return geo;
}

Geometry& normalizeNormals(Geometry& geo)
{
	for (auto& n : geo.normals)
		n = normalize(n);
	return geo;
}

Geometry& flipNormals(Geometry& geo)
{
	// TODO: Should probably flip face winding too
	for (auto& n : geo.normals)
		n = -n;
	return geo;
}

Geometry& generateBoxTexcoords(Geometry& geo)
{
	if (geo.positions.empty())
		return geo;
	if (geo.normals.empty())
		calculateNormals(geo);
	geo.texcoords.resize(geo.positions.size());
	Bounds bb = bounds(geo);
	vec3 size = bb.size();
	for (int i = 0; i < geo.positions.size(); ++i) {
		vec3 p = geo.positions[i] - bb.min;
		vec3 n = geo.normals[i];
		vec2& uv = geo.texcoords[i];
		// Determine dominant cardinal direction
		vec3 an = abs(n);
		if (an.x >= an.y && an.x >= an.z) { // facing x
			uv = { p.z / size.z, p.y / size.y };
			if (n.x > 0) uv.x = 1.f - uv.x;
		} else if (an.y > an.z) { // facing y
			uv = { p.x / size.x, 1.0f - p.z / size.z };
			if (n.y < 0) uv.x = 1.f - uv.x;
		} else { // facing z
			uv = { p.x / size.x, p.y / size.y };
			if (n.z < 0) uv.x = 1.f - uv.x;
		}
		uv.y = FlipYForDX(uv.y, UvStyle);
	}
	return geo;
}


bool writeObj(const Geometry& geo, const char* filename, bool optimize)
{
	std::ofstream f;
	f.open(filename, std::fstream::binary | std::fstream::out | std::fstream::trunc); // We don't want line ending handling...
	if (!f.good())
		return false;

	f << "# Generated by GenGeo" << std::endl;

	char faceStyle = 0;
	if (!geo.texcoords.empty())
		faceStyle |= 1;
	if (!geo.normals.empty())
		faceStyle |= 2;

	if (!geo.colors.empty())
		optimize = false; // TODO: Optimize not currently supported for vertex colors

	if (optimize) {
		// Spend some effort to only output the minimal amount of data
		constexpr double tolerance = 0.00001;
		auto quant3 = [=](vec3 v) {
			v.x = (float)(std::round((double)v.x / tolerance) * tolerance);
			v.y = (float)(std::round((double)v.y / tolerance) * tolerance);
			v.z = (float)(std::round((double)v.z / tolerance) * tolerance);
			return v;
		};
		auto quant2 = [=](vec2 v) {
			v.x = (float)(std::round((double)v.x / tolerance) * tolerance);
			v.y = (float)(std::round((double)v.y / tolerance) * tolerance);
			return v;
		};

		std::unordered_map<vec3, uint> hashmapPos;
		int index = 1;
		for (vec3 v : geo.positions) {
			v = quant3(v);
			uint newIndex = hashmapPos.try_emplace(v, index).first->second;
			if (index == newIndex) {
				f << tfm::format("v %f %f %f\n", v.x, v.y, v.z);
				index++;
			}
		}

		std::unordered_map<vec2, uint> hashmapUv;
		index = 1;
		for (vec2 uv : geo.texcoords) {
			uv = quant2(uv);
			uint newIndex = hashmapUv.try_emplace(uv, index).first->second;
			if (index == newIndex) {
				f << tfm::format("vt %f %f\n", uv.x, uv.y);
				index++;
			}
		}

		std::unordered_map<vec3, uint> hashmapNorm;
		index = 1;
		for (vec3 n : geo.normals) {
			n = quant3(n);
			uint newIndex = hashmapNorm.try_emplace(n, index).first->second;
			if (index == newIndex) {
				f << tfm::format("vn %f %f %f\n", n.x, n.y, n.z);
				index++;
			}
		}

		auto ip = [&](int index) {
			vec3 v = quant3(geo.positions[index]);
			return hashmapPos[v];
		};
		auto it = [&](int index) {
			vec2 v = quant2(geo.texcoords[index]);
			return hashmapUv[v];
		};
		auto in = [&](int index) {
			vec3 v = quant3(geo.normals[index]);
			return hashmapNorm[v];
		};

		if (geo.triangles.empty()) {
			for (const auto& face : geo.quads) {
				switch (faceStyle) {
					case 0: f << tfm::format("f %d %d %d %d\n", ip(face.a), ip(face.b), ip(face.c), ip(face.d)); break;
					case 1: f << tfm::format("f %d/%d %d/%d %d/%d %d/%d\n", ip(face.a), it(face.a), ip(face.b), it(face.b), ip(face.c), it(face.c), ip(face.d), it(face.d)); break;
					case 2: f << tfm::format("f %d//%d %d//%d %d//%d %d//%d\n", ip(face.a), in(face.a), ip(face.b), in(face.b), ip(face.c), in(face.c), ip(face.d), in(face.d)); break;
					case 3: f << tfm::format("f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d\n",
						ip(face.a), it(face.a), in(face.a),
						ip(face.b), it(face.b), in(face.b),
						ip(face.c), it(face.c), in(face.c),
						ip(face.d), it(face.d), in(face.d)); break;
				}
			}
		} else {
			for (const auto& face : geo.triangles) {
				switch (faceStyle) {
					case 0: f << tfm::format("f %d %d %d\n", ip(face.a), ip(face.b), ip(face.c)); break;
					case 1: f << tfm::format("f %d/%d %d/%d %d/%d\n", ip(face.a), it(face.a), ip(face.b), it(face.b), ip(face.c), it(face.c)); break;
					case 2: f << tfm::format("f %d//%d %d//%d %d//%d\n", ip(face.a), in(face.a), ip(face.b), in(face.b), ip(face.c), in(face.c)); break;
					case 3: f << tfm::format("f %d/%d/%d %d/%d/%d %d/%d/%d\n",
						ip(face.a), it(face.a), in(face.a),
						ip(face.b), it(face.b), in(face.b),
						ip(face.c), it(face.c), in(face.c)); break;
				}
			}
		}
	// Just dump everything as is to the file
	} else {
		// Non standard, but somewhat supported extension for vertex colors
		// http://paulbourke.net/dataformats/obj/colour.html
		// Alpha omitted for better(?) compatibility
		if (geo.colors.size() == geo.positions.size()) {
			for (uint i = 0; i < geo.positions.size(); ++i) {
				vec3 v = geo.positions[i];
				Color c = geo.colors[i];
				f << tfm::format("v %f %f %f %f %f %f\n", v.x, v.y, v.z, c.r, c.g, c.b);
			}
		} else {
			for (const auto& v : geo.positions)
				f << tfm::format("v %f %f %f\n", v.x, v.y, v.z);
		}

		for (const auto& uv : geo.texcoords)
			f << tfm::format("vt %f %f\n", uv.x, uv.y);

		for (const auto& n : geo.normals)
			f << tfm::format("vn %f %f %f\n", n.x, n.y, n.z);

		if (geo.triangles.empty()) {
			for (const auto& face : geo.quads) {
				int a = face.a + 1, b = face.b + 1, c = face.c + 1, d = face.d + 1;
				switch (faceStyle) {
					case 0: f << tfm::format("f %d %d %d %d\n", a, b, c, d); break;
					case 1: f << tfm::format("f %d/%d %d/%d %d/%d %d/%d\n", a, a, b, b, c, c, d, d); break;
					case 2: f << tfm::format("f %d//%d %d//%d %d//%d %d//%d\n", a, a, b, b, c, c, d, d); break;
					case 3: f << tfm::format("f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d\n", a, a, a, b, b, b, c, c, c, d, d, d); break;
				}
			}
		} else {
			for (const auto& face : geo.triangles) {
				int a = face.a + 1, b = face.b + 1, c = face.c + 1;
				switch (faceStyle) {
					case 0: f << tfm::format("f %d %d %d\n", a, b, c); break;
					case 1: f << tfm::format("f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c); break;
					case 2: f << tfm::format("f %d//%d %d//%d %d//%d\n", a, a, b, b, c, c); break;
					case 3: f << tfm::format("f %d/%d/%d %d/%d/%d %d/%d/%d\n", a, a, a, b, b, b, c, c, c); break;
				}
			}
		}
	}

	return true;
}

bool writeStl(const Geometry& inGeo, const char* filename, MeshFormatStyle formatStyle)
{
	std::ofstream f;
	f.open(filename, std::fstream::out | std::fstream::trunc | std::ios::binary);
	if (!f.good())
		return false;

	Geometry geo = inGeo;
	triangulate(geo); // Only triangles supported
	rotatex(geo, 90); // TODO: Allow configuration (but convention seems to be z-up for STL / 3d printing)
	if (geo.normals.empty())
		calculateNormals(geo);

	if (formatStyle == MeshFormatStyle::Binary) {
		const char header[80] = {0};
		f.write(header, sizeof(header));
		unsigned int numTris = geo.triangles.size();
		f.write((const char*)&numTris, 4); // UINT32 - Triangle count
		for (auto tri : geo.triangles) {
			vec3 n = normalize(geo.normals[tri.a] + geo.normals[tri.b] + geo.normals[tri.c]);
			vec3 p1 = geo.positions[tri.a];
			vec3 p2 = geo.positions[tri.b];
			vec3 p3 = geo.positions[tri.c];
			f.write((const char*)&n, sizeof(vec3));
			f.write((const char*)&p1, sizeof(vec3));
			f.write((const char*)&p2, sizeof(vec3));
			f.write((const char*)&p3, sizeof(vec3));
			f.put(0); // UINT16 – Attribute byte count
			f.put(0);
		}
	} else {
		f << "solid \n"; // Could append some kind of name here, note space required
		for (auto tri : geo.triangles) {
			vec3 n = normalize(geo.normals[tri.a] + geo.normals[tri.b] + geo.normals[tri.c]);
			vec3 p1 = geo.positions[tri.a];
			vec3 p2 = geo.positions[tri.b];
			vec3 p3 = geo.positions[tri.c];
			f << "facet normal " << n.x << " " << n.y << " " << n.z << "\n";
			f << "\touter loop\n";
			f << "\t\t vertex " << p1.x << " " << p1.y << " " << p1.z << "\n";
			f << "\t\t vertex " << p2.x << " " << p2.y << " " << p2.z << "\n";
			f << "\t\t vertex " << p3.x << " " << p3.y << " " << p3.z << "\n";
			f << "\tendloop\n";
			f << "endfacet\n";
		}
		f << "endsolid \n"; // Name also here
		f << std::flush;
	}

	return true;
}


namespace sdf {

template<typename T>
inline void push(std::vector<T>& stack, T value) {
	stack.push_back(value);
}

template<typename T>
inline void push_front(std::vector<T>& stack, T value) {
	stack.insert(stack.begin(), value);
}

inline void push(std::vector<float>& stack, vec3 v) {
	stack.push_back(v.x);
	stack.push_back(v.y);
	stack.push_back(v.z);
}

inline float popFloat(std::vector<float>& stack) {
	float value = stack.back();
	stack.pop_back();
	return value;
}

// If used check component popping order
/*inline vec3 popVec3(std::vector<float>& stack) {
	vec3 ret;
	ret.x = popFloat(stack);
	ret.y = popFloat(stack);
	ret.z = popFloat(stack);
	return ret;
}*/

inline float getFloat(float*& paramsPtr) {
	float ret = *paramsPtr;
	++paramsPtr;
	return ret;
}

inline vec3 getVec3(float*& paramsPtr) {
	vec3 ret = { paramsPtr[0], paramsPtr[1], paramsPtr[2] };
	paramsPtr += 3;
	return ret;
}

template<typename T, uint N>
struct StaticStack {
	static constexpr uint MAX_STACK = N;
	int stackPtr = -1;
	T stack[MAX_STACK];

	void push(T value) { stack[++stackPtr] = value; }
	T pop() { return stack[stackPtr--]; }
	bool empty() const { return stackPtr < 0; }
};

float SDF::operator()(vec3 pos)
{
	StaticStack<float, 128> stack;
	StaticStack<vec3, 16> posStack;
	int count = code.size();
	float* paramsPtr = &params.front();
	int geoIndex = 0;
	for (int i = 0; i < count; ++i) {
		auto op = code[i];
		switch (op) {
			case NOOP:
				break;
			case PRIM_SPHERE: {
				const vec3 p = translate(pos, getVec3(paramsPtr));
				stack.push(sdf::sphere(p, getFloat(paramsPtr)));
				break;
			}
			case PRIM_BOX: {
				const vec3 p = translate(pos, getVec3(paramsPtr));
				stack.push(sdf::box(p, getVec3(paramsPtr)));
				break;
			}
			case PRIM_PLANE: {
				const vec3 p = translate(pos, getVec3(paramsPtr));
				const vec3 n = getVec3(paramsPtr);
				stack.push(sdf::plane(p, n, getFloat(paramsPtr)));
				break;
			}
			case PRIM_CAPSULE: {
				const vec3 a = getVec3(paramsPtr);
				const vec3 b = getVec3(paramsPtr);
				const float r = getFloat(paramsPtr);
				stack.push(sdf::capsule(pos, a, b, r));
				break;
			}
			case PRIM_CAPSULEX: {
				const vec3 p = translate(pos, getVec3(paramsPtr));
				const float h = getFloat(paramsPtr);
				const float r = getFloat(paramsPtr);
				stack.push(sdf::capsulex(p, h, r));
				break;
			}
			case PRIM_CAPSULEY: {
				const vec3 p = translate(pos, getVec3(paramsPtr));
				const float h = getFloat(paramsPtr);
				const float r = getFloat(paramsPtr);
				stack.push(sdf::capsuley(p, h, r));
				break;
			}
			case PRIM_CAPSULEZ: {
				const vec3 p = translate(pos, getVec3(paramsPtr));
				const float h = getFloat(paramsPtr);
				const float r = getFloat(paramsPtr);
				stack.push(sdf::capsulez(p, h, r));
				break;
			}
			case PRIM_CYLINDER: {
				const vec3 a = getVec3(paramsPtr);
				const vec3 b = getVec3(paramsPtr);
				const float r = getFloat(paramsPtr);
				stack.push(sdf::cylinder(pos, a, b, r));
				break;
			}
			case PRIM_CYLINDERX: {
				const vec3 p = translate(pos, getVec3(paramsPtr));
				const float h = getFloat(paramsPtr);
				const float r = getFloat(paramsPtr);
				stack.push(sdf::cylinderx(p, h, r));
				break;
			}
			case PRIM_CYLINDERY: {
				const vec3 p = translate(pos, getVec3(paramsPtr));
				const float h = getFloat(paramsPtr);
				const float r = getFloat(paramsPtr);
				stack.push(sdf::cylindery(p, h, r));
				break;
			}
			case PRIM_CYLINDERZ: {
				const vec3 p = translate(pos, getVec3(paramsPtr));
				const float h = getFloat(paramsPtr);
				const float r = getFloat(paramsPtr);
				stack.push(sdf::cylinderz(p, h, r));
				break;
			}
			case OP_UNION: {
				float b = stack.pop();
				float a = stack.pop();
				stack.push(sdf::opUnion(a, b));
				break;
			}
			case OP_SUBTRACT: {
				float b = stack.pop();
				float a = stack.pop();
				stack.push(sdf::opSubtract(a, b));
				break;
			}
			case OP_INTERSECT: {
				float b = stack.pop();
				float a = stack.pop();
				stack.push(sdf::opIntersect(a, b));
				break;
			}
			case OP_UNION_SMOOTH: {
				float b = stack.pop();
				float a = stack.pop();
				stack.push(sdf::opSmoothUnion(a, b, getFloat(paramsPtr)));
				break;
			}
			case OP_SUBTRACT_SMOOTH: {
				float b = stack.pop();
				float a = stack.pop();
				stack.push(sdf::opSmoothSubtraction(a, b, getFloat(paramsPtr)));
				break;
			}
			case OP_INTERSECT_SMOOTH: {
				float b = stack.pop();
				float a = stack.pop();
				stack.push(sdf::opSmoothIntersection(a, b, getFloat(paramsPtr)));
				break;
			}
			case OP_ROUND: {
				stack.push(sdf::opRound(stack.pop(), getFloat(paramsPtr)));
				break;
			}
			case OP_ONION: {
				stack.push(sdf::opOnion(stack.pop(), getFloat(paramsPtr)));
				break;
			}
			case OP_SYMMETRY: {
				if (getFloat(paramsPtr) > 0.f)
					pos.x = abs(pos.x);
				if (getFloat(paramsPtr) > 0.f)
					pos.y = abs(pos.y);
				if (getFloat(paramsPtr) > 0.f)
					pos.z = abs(pos.z);
				break;
			}

			case PUSH_POS: {
				posStack.push(pos);
				break;
			}
			case POP_POS: {
				pos = posStack.pop();
				break;
			}
		}
	}
	assert(stack.stackPtr == 0);
	assert(paramsPtr == &params.back() + 1);
	return stack.pop();
}

SDF& SDF::combine(const SDF& other, SDFInstruction op)
{
	if (other.empty())
		return *this;
	bool meEmpty = empty();
	code.insert(code.end(), other.code.begin(), other.code.end());
	if (!meEmpty) // Op is ignored if we are combining with empty
		code.push_back(op);
	params.insert(params.end(), other.params.begin(), other.params.end());
	geos.insert(geos.end(), other.geos.begin(), other.geos.end());
	geoBounds.insert(geoBounds.end(), other.geoBounds.begin(), other.geoBounds.end());
	return *this;
}

SDF& SDF::round(float r)
{
	code.push_back(OP_ROUND);
	push(params, r);
	return *this;
}

SDF& SDF::onion(float thickness)
{
	code.push_back(OP_ONION);
	push(params, thickness);
	return *this;
}

SDF& SDF::symmetry(bool x, bool y, bool z)
{
	code.insert(code.begin(), OP_SYMMETRY); // Symmetry needs to affect sampling pos so go to beginning to adjust it
	push_front(params, z ? 1.f : 0.f); // Parameters will also need to be at the beginning
	push_front(params, y ? 1.f : 0.f); // Reverse order for parameters due to push_front
	push_front(params, x ? 1.f : 0.f);
	code.insert(code.begin(), PUSH_POS); // We want to be able to restore pos, so save it as the very first thing
	code.push_back(POP_POS); // At last, restore pos for everything that comes after
	return *this;
}

SDF& SDF::sphere(vec3 pos, float r)
{
	code.push_back(PRIM_SPHERE);
	push(params, pos);
	push(params, r);
	return *this;
}

SDF& SDF::box(vec3 pos, vec3 bounds)
{
	code.push_back(PRIM_BOX);
	push(params, pos);
	push(params, bounds);
	return *this;
}

SDF& SDF::plane(vec3 pos, vec3 n, float d)
{
	code.push_back(PRIM_PLANE);
	push(params, pos);
	push(params, n);
	push(params, d);
	return *this;
}

SDF& SDF::capsule(vec3 a, vec3 b, float r)
{
	code.push_back(PRIM_CAPSULE);
	push(params, a);
	push(params, b);
	push(params, r);
	return *this;
}

SDF& SDF::capsulex(vec3 pos, float h, float r)
{
	code.push_back(PRIM_CAPSULEX);
	push(params, pos);
	push(params, h);
	push(params, r);
	return *this;
}

SDF& SDF::capsuley(vec3 pos, float h, float r)
{
	code.push_back(PRIM_CAPSULEY);
	push(params, pos);
	push(params, h);
	push(params, r);
	return *this;
}

SDF& SDF::capsulez(vec3 pos, float h, float r)
{
	code.push_back(PRIM_CAPSULEZ);
	push(params, pos);
	push(params, h);
	push(params, r);
	return *this;
}

SDF& SDF::cylinder(vec3 a, vec3 b, float r)
{
	code.push_back(PRIM_CYLINDER);
	push(params, a);
	push(params, b);
	push(params, r);
	return *this;
}

SDF& SDF::cylinderx(vec3 pos, float h, float r)
{
	code.push_back(PRIM_CYLINDERX);
	push(params, pos);
	push(params, h);
	push(params, r);
	return *this;
}

SDF& SDF::cylindery(vec3 pos, float h, float r)
{
	code.push_back(PRIM_CYLINDERY);
	push(params, pos);
	push(params, h);
	push(params, r);
	return *this;
}

SDF& SDF::cylinderz(vec3 pos, float h, float r)
{
	code.push_back(PRIM_CYLINDERZ);
	push(params, pos);
	push(params, h);
	push(params, r);
	return *this;
}

std::vector<float> slice(vec2 size, vec3 topLeft, vec3 topRight, vec3 bottomLeft, vec3 bottomRight, float offset, SDFFunction sdf)
{
	std::vector<float> slice;
	slice.resize(size.x * size.y);
	int w = size.x, h = size.y;
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			float u = x / (size.x - 1); // Don't use w and h because they are int
			float v = y / (size.y - 1);
			vec3 top = lerp(topLeft, topRight, u);
			vec3 bottom = lerp(bottomLeft, bottomRight, u);
			vec3 pos = lerp(top, bottom, v);
			slice[y * w + x] = sdf(pos) + offset;
		}
	}
	return slice;
}

} // namespace sdf


VoxelGrid& box(VoxelGrid& voxelGrid, vec3 start, vec3 end, VoxelGrid::cell_t cell)
{
	voxelGrid.visit([cell](int x, int y, int z) {
		return cell;
	});
	return voxelGrid;
}

Geometry polygonize(const VoxelGrid& voxelGrid, float voxelSize)
{
	Geometry geo;
	vec3 offset = vec3(voxelGrid.xsize, voxelGrid.ysize, voxelGrid.zsize) * -0.5f - 0.5f;
	voxelGrid.visit([&voxelGrid, &geo, offset, voxelSize](int x, int y, int z, VoxelGrid::cell_t cell) {
		if (!cell)
			return;
		auto prevcellx = voxelGrid.getSafe(x - 1, y, z);
		auto prevcelly = voxelGrid.getSafe(x, y - 1, z);
		auto prevcellz = voxelGrid.getSafe(x, y, z - 1);
		auto nextcellx = voxelGrid.getSafe(x + 1, y, z);
		auto nextcelly = voxelGrid.getSafe(x, y + 1, z);
		auto nextcellz = voxelGrid.getSafe(x, y, z + 1);
		auto handleSide = [=, &geo](Geometry plane, vec3 sideOffset, vec3 rot) {
			// TODO: This could be optimized...
			scale(plane, vec3(voxelSize));
			rotate(plane, rot);
			translate(plane, (vec3(x, y, z) + offset + sideOffset * 0.5f) * voxelSize);
			concat(geo, plane);
		};
		if (prevcellx == 0) handleSide(planeyz(), { -1,  0,  0 }, { 0, 180, 0 });
		if (prevcelly == 0) handleSide(planexz(), { 0, -1,  0 }, { 0, 0, 180 });
		if (prevcellz == 0) handleSide(planexy(), { 0,  0, -1 }, { 0, 180, 0 });
		if (nextcellx == 0) handleSide(planeyz(), { 1,  0,  0 }, {});
		if (nextcelly == 0) handleSide(planexz(), { 0,  1,  0 }, {});
		if (nextcellz == 0) handleSide(planexy(), { 0,  0,  1 }, {});
	});
	return geo;
}

// Linearly interpolate the position where an isosurface cuts
// an edge between two vertices, each with their own scalar value
static constexpr vec3 vertexInterp(float isolevel, vec3 p1, vec3 p2, float valp1, float valp2)
{
	if (abs(isolevel - valp1) < 0.0001f)
		return p1;
	if (abs(isolevel - valp2) < 0.0001f)
		return p2;
	if (abs(valp1 - valp2) < 0.0001f)
		return p1;
	float mu = (isolevel - valp1) / (valp2 - valp1);
	return {
		p1.x + mu * (p2.x - p1.x),
		p1.y + mu * (p2.y - p1.y),
		p1.z + mu * (p2.z - p1.z)
	};
}

namespace sdf {

Geometry polygonize(vec3 start, vec3 end, vec3 step, SDFFunction sdf)
{
	static constexpr vec3 offsets[8] = {
		{-1,  1, -1},
		{ 1,  1, -1},
		{ 1, -1, -1},
		{-1, -1, -1},
		{-1,  1,  1},
		{ 1,  1,  1},
		{ 1, -1,  1},
		{-1, -1,  1}
	};
	Geometry geo;
	if (step.x == 0)
		step.x = 1000000.f;
	if (step.y == 0)
		step.y = 1000000.f;
	if (step.z == 0)
		step.z = 1000000.f;
	const vec3 offsetMult = step * 0.5f;
	const float isolevel = 0;
	for (float z = start.z; z < end.z; z += step.z) {
		for (float y = start.y; y < end.y; y += step.y) {
			for (float x = start.x; x < end.x; x += step.x) {
				// Determine the index into the edge table which
				// tells us which vertices are inside of the surface
				int cubeindex = 0;
				float values[8];
				vec3 points[8];
				const vec3 p = vec3(x, y, z);
				for (int i = 0; i < 8; i++) {
					points[i] = p + offsets[i] * offsetMult;
					values[i] = sdf(points[i]);
					if (values[i] < isolevel)
						cubeindex |= (1 << i);
				}

				// Cube is entirely in/out of the surface
				if (EdgeTable[cubeindex] == 0)
					continue;

				// Find the vertices where the surface intersects the cube
				vec3 vertlist[12];
				if (EdgeTable[cubeindex] & 1)
					vertlist[0] = vertexInterp(isolevel, points[0], points[1], values[0], values[1]);
				if (EdgeTable[cubeindex] & 2)
					vertlist[1] = vertexInterp(isolevel, points[1], points[2], values[1], values[2]);
				if (EdgeTable[cubeindex] & 4)
					vertlist[2] = vertexInterp(isolevel, points[2], points[3], values[2], values[3]);
				if (EdgeTable[cubeindex] & 8)
					vertlist[3] = vertexInterp(isolevel, points[3], points[0], values[3], values[0]);
				if (EdgeTable[cubeindex] & 16)
					vertlist[4] = vertexInterp(isolevel, points[4], points[5], values[4], values[5]);
				if (EdgeTable[cubeindex] & 32)
					vertlist[5] = vertexInterp(isolevel, points[5], points[6], values[5], values[6]);
				if (EdgeTable[cubeindex] & 64)
					vertlist[6] = vertexInterp(isolevel, points[6], points[7], values[6], values[7]);
				if (EdgeTable[cubeindex] & 128)
					vertlist[7] = vertexInterp(isolevel, points[7], points[4], values[7], values[4]);
				if (EdgeTable[cubeindex] & 256)
					vertlist[8] = vertexInterp(isolevel, points[0], points[4], values[0], values[4]);
				if (EdgeTable[cubeindex] & 512)
					vertlist[9] = vertexInterp(isolevel, points[1], points[5], values[1], values[5]);
				if (EdgeTable[cubeindex] & 1024)
					vertlist[10] = vertexInterp(isolevel, points[2], points[6], values[2], values[6]);
				if (EdgeTable[cubeindex] & 2048)
					vertlist[11] = vertexInterp(isolevel, points[3], points[7], values[3], values[7]);

				// Create the triangle
				for (int i = 0; TriTable[cubeindex][i] != -1; i += 3) {
					const vec3 a = vertlist[TriTable[cubeindex][i  ]];
					const vec3 b = vertlist[TriTable[cubeindex][i+1]];
					const vec3 c = vertlist[TriTable[cubeindex][i+2]];
					geo.positions.push_back(a);
					geo.positions.push_back(b);
					geo.positions.push_back(c);
					const vec3 n = normalize(cross(b - a, c - a));
					geo.normals.push_back(n);
					geo.normals.push_back(n);
					geo.normals.push_back(n);
					const int s = geo.positions.size();
					geo.triangles.emplace_back(s - 3, s - 2, s - 1);
				}
			}
		}
	}
	return geo;
}

} // namespace sdf


Geometry polygonizeSmooth(const VoxelGrid& voxelGrid, float isolevel)
{
	static constexpr vec3 offsets[8] = {
		{-1,  1, -1},
		{ 1,  1, -1},
		{ 1, -1, -1},
		{-1, -1, -1},
		{-1,  1,  1},
		{ 1,  1,  1},
		{ 1, -1,  1},
		{-1, -1,  1}
	};
	Geometry geo;
	voxelGrid.visit([&voxelGrid, &geo, isolevel](int x, int y, int z, VoxelGrid::cell_t cell) {
		// Determine the index into the edge table which
		// tells us which vertices are inside of the surface
		int cubeindex = 0;
		float values[8];
		vec3 points[8];
		const vec3 gridCenter = vec3(voxelGrid.xsize, voxelGrid.ysize, voxelGrid.zsize) * 0.5f;
		const vec3 p = vec3(x, y, z);
		//auto current = voxelGrid.getSafe(p);
		for (int i = 0; i < 8; i++) {
			points[i] = p + offsets[i] * 0.5f;
			//auto corner = voxelGrid.getSafe(p + offsets[i]);
			// TODO: This is actually just sphere, not using voxel info
			values[i] = length(gridCenter - points[i]) - (gridCenter.x * 0.9f);
			if (values[i] < isolevel)
				cubeindex |= (1 << i);
		}

		// Cube is entirely in/out of the surface
		if (EdgeTable[cubeindex] == 0)
			return;

		// Find the vertices where the surface intersects the cube
		vec3 vertlist[12];
		if (EdgeTable[cubeindex] & 1)
			vertlist[0] = vertexInterp(isolevel, points[0], points[1], values[0], values[1]);
		if (EdgeTable[cubeindex] & 2)
			vertlist[1] = vertexInterp(isolevel, points[1], points[2], values[1], values[2]);
		if (EdgeTable[cubeindex] & 4)
			vertlist[2] = vertexInterp(isolevel, points[2], points[3], values[2], values[3]);
		if (EdgeTable[cubeindex] & 8)
			vertlist[3] = vertexInterp(isolevel, points[3], points[0], values[3], values[0]);
		if (EdgeTable[cubeindex] & 16)
			vertlist[4] = vertexInterp(isolevel, points[4], points[5], values[4], values[5]);
		if (EdgeTable[cubeindex] & 32)
			vertlist[5] = vertexInterp(isolevel, points[5], points[6], values[5], values[6]);
		if (EdgeTable[cubeindex] & 64)
			vertlist[6] = vertexInterp(isolevel, points[6], points[7], values[6], values[7]);
		if (EdgeTable[cubeindex] & 128)
			vertlist[7] = vertexInterp(isolevel, points[7], points[4], values[7], values[4]);
		if (EdgeTable[cubeindex] & 256)
			vertlist[8] = vertexInterp(isolevel, points[0], points[4], values[0], values[4]);
		if (EdgeTable[cubeindex] & 512)
			vertlist[9] = vertexInterp(isolevel, points[1], points[5], values[1], values[5]);
		if (EdgeTable[cubeindex] & 1024)
			vertlist[10] = vertexInterp(isolevel, points[2], points[6], values[2], values[6]);
		if (EdgeTable[cubeindex] & 2048)
			vertlist[11] = vertexInterp(isolevel, points[3], points[7], values[3], values[7]);

		// Create the triangle
		for (int i = 0; TriTable[cubeindex][i] != -1; i += 3) {
			const vec3 a = vertlist[TriTable[cubeindex][i  ]];
			const vec3 b = vertlist[TriTable[cubeindex][i+1]];
			const vec3 c = vertlist[TriTable[cubeindex][i+2]];
			geo.positions.push_back(a);
			geo.positions.push_back(b);
			geo.positions.push_back(c);
			const vec3 n = normalize(cross(b - a, c - a));
			geo.normals.push_back(n);
			geo.normals.push_back(n);
			geo.normals.push_back(n);
			const int s = geo.positions.size();
			geo.triangles.emplace_back(s - 3, s - 2, s - 1);
		}
	});
	return geo;
}

} // namespace gengeo
