#include "gengeo.hpp"
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <ctime>
#include <random>
#include <unordered_map>

namespace gengeo {

#include "marchingcubes.h"

static constexpr float PI = 3.1415926535f;
static constexpr float TWOPI = 2 * PI;
//static constexpr float PI2 = PI * 0.5f;
//static constexpr float PI4 = PI * 0.25f;
static constexpr float DEG_TO_RAD = PI / 180.f;
//static constexpr float RAD_TO_DEG = 180.f / PI;

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
		geo.texcoords.emplace_back(0.f, 0.f);
		geo.texcoords.emplace_back(1.f, 0.f);
		geo.texcoords.emplace_back(1.f, 1.f);
		geo.texcoords.emplace_back(0.f, 1.f);
		geo.normals.emplace_back(normals[i]);
		geo.normals.emplace_back(normals[i]);
		geo.normals.emplace_back(normals[i]);
		geo.normals.emplace_back(normals[i]);
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
	geo.texcoords.emplace_back(0.f, 0.f);
	geo.texcoords.emplace_back(1.f, 0.f);
	geo.texcoords.emplace_back(1.f, 1.f);
	geo.texcoords.emplace_back(0.f, 1.f);
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
	geo.texcoords.emplace_back(0.f, 0.f);
	geo.texcoords.emplace_back(1.f, 0.f);
	geo.texcoords.emplace_back(1.f, 1.f);
	geo.texcoords.emplace_back(0.f, 1.f);
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
	geo.texcoords.emplace_back(0.f, 0.f);
	geo.texcoords.emplace_back(1.f, 0.f);
	geo.texcoords.emplace_back(1.f, 1.f);
	geo.texcoords.emplace_back(0.f, 1.f);
	geo.normals.assign(4, { 1.f, 0.f, 0.f });
	return geo;
}

Geometry heightmap(float* heights, int w, int h, vec3 scale)
{
	int numVerts = w * h;
	int lasti = w - 1;
	int lastj = h - 1;
	Geometry geo;
	geo.positions.resize(numVerts);
	geo.texcoords.resize(numVerts);
	geo.normals.resize(numVerts);
	geo.quads.reserve(numVerts);
	for (int j = 0; j < h; ++j) {
		for (int i = 0; i < w; ++i) {
			// Create the vertex
			float x = i - w * 0.5f;
			float z = j - h * 0.5f;
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

void bounds(const Geometry& geo, vec3& lower, vec3& upper)
{
	vec3 lbound(std::numeric_limits<float>::max());
	vec3 ubound(-std::numeric_limits<float>::max());
	for (auto& pos : geo.positions) {
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
	// TODO: Normals (non-uniform, mirroring)
	for (auto& pos : geo.positions)
		pos *= s;
	return geo;
}

Geometry& translate(Geometry& geo, vec3 v)
{
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

Geometry& concat(Geometry& geo, const Geometry& other)
{
	uint indexOffset = geo.positions.size();
	geo.positions.insert(geo.positions.end(), other.positions.begin(), other.positions.end());
	geo.texcoords.insert(geo.texcoords.end(), other.texcoords.begin(), other.texcoords.end());
	geo.normals.insert(geo.normals.end(), other.normals.begin(), other.normals.end());

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
		subdivided.quads.reserve(geo.quads.size() * 4);
		for (const auto& quad : geo.quads) {
			uint base = subdivided.positions.size();
			subdiv(quad, subdivided.positions, geo.positions);
			subdiv(quad, subdivided.texcoords, geo.texcoords);
			subdiv(quad, subdivided.normals, geo.normals);
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
}

Geometry& weld(Geometry& geo, float maxNormalAngleDeg, double tolerance)
{
	// Round to tolerance
	for (auto& pos : geo.positions) {
		pos.x = (float)(std::round((double)pos.x / tolerance) * tolerance);
		pos.y = (float)(std::round((double)pos.y / tolerance) * tolerance);
		pos.z = (float)(std::round((double)pos.z / tolerance) * tolerance);
	}
	float cosTheta = std::cos(maxNormalAngleDeg * DEG_TO_RAD);
	std::unordered_map<vec3, uint> hashmap;
	Geometry welded;
	welded.positions.reserve(geo.positions.size());
	welded.texcoords.reserve(geo.texcoords.size());
	welded.normals.reserve(geo.normals.size());

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
	weld(geo, 180);  // Need to weld everything to avoid cracks...
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
	if (geo.texcoords.empty())
		return geo;
	weld(geo, 180);  // Need to weld everything to avoid cracks...
	if (geo.normals.empty())
		calculateNormals(geo);
	auto& pos = geo.positions;
	auto& uv = geo.texcoords;
	auto& n = geo.normals;
	for (uint i = 0; i < pos.size(); i++)
		pos[i] += n[i] * sampleTexture(uv[i], displacementTexture, w, h) * height;
	return calculateNormals(geo);
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
			vec3 normal = triangleNormal(positions[quad.a], positions[quad.b], positions[quad.c]);
			normals[quad.a] += normal;
			normals[quad.b] += normal;
			normals[quad.c] += normal;
			normals[quad.d] += normal;
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

bool writeObj(const Geometry& geo, const char* filename)
{
	std::ofstream f;
	f.open(filename, std::fstream::out | std::fstream::trunc);
	if (!f.good())
		return false;

	f << "# Generated by GenGeo" << std::endl;

	for (const auto& v : geo.positions)
		f << "v " << v.x << " " << v.y << " " << v.z << std::endl;

	for (const auto& uv : geo.texcoords)
		f << "vt " << uv.x << " " << uv.y << std::endl;

	for (const auto& n : geo.normals)
		f << "vn " << n.x << " " << n.y << " " << n.z << std::endl;

	// TODO: Handle missing normals/texcoords

	if (geo.triangles.empty()) {
		for (const auto& face : geo.quads) {
			int a = face.a + 1;
			int b = face.b + 1;
			int c = face.c + 1;
			int d = face.d + 1;
			f << "f " << a << "/" << a << "/" << a
				<< " " << b << "/" << b << "/" << b
				<< " " << c << "/" << c << "/" << c
				<< " " << d << "/" << d << "/" << d << std::endl;
		}
	} else {
		for (const auto& face : geo.triangles) {
			int a = face.a + 1;
			int b = face.b + 1;
			int c = face.c + 1;
			f << "f " << a << "/" << a << "/" << a
				<< " " << b << "/" << b << "/" << b
				<< " " << c << "/" << c << "/" << c << std::endl;
		}
	}

	return true;
}


VoxelGrid& box(VoxelGrid& voxelGrid, vec3 start, vec3 end, VoxelGrid::cell_t cell)
{
	voxelGrid.visit([cell](int x, int y, int z) {
		return cell;
	});
	return voxelGrid;
}

// TODO: This is super quick and dirty, no internal face removal
Geometry polygonize(const VoxelGrid& voxelGrid)
{
	Geometry geo;
	Geometry voxel = cube();
	voxelGrid.visit([&](int x, int y, int z, VoxelGrid::cell_t cell) {
		if (!cell)
			return;
		translate(center(voxel), vec3(x, y, z));
		concat(geo, voxel);
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
					points[i] = p + offsets[i] * 0.5f;
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
			// TODO: This is actually just spehere, not using voxel info
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
