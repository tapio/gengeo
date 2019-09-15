#include "gengeo.hpp"
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>

namespace gengeo {

static constexpr float PI = 3.1415926535f;
static constexpr float TWOPI = 2 * PI;
static constexpr float PI2 = PI * 0.5f;
static constexpr float PI4 = PI * 0.25f;

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
		geo.texcoords.emplace_back(0, 0);
		geo.texcoords.emplace_back(1, 0);
		geo.texcoords.emplace_back(1, 1);
		geo.texcoords.emplace_back(0, 1);
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
	subdivide(scale(geo, 2), subdivs);
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
	geo.positions.emplace_back(0, 0, 0);
	geo.texcoords.emplace_back(0.5f, 0.5f);
	float angleStep = TWOPI / points;
	for (int i = 0; i < points; i++) {
		geo.triangles.emplace_back(0, ((i+1) % points) + 1, i+1);
		float x = std::cos(i * angleStep);
		float y = std::sin(i * angleStep);
		float u = (x + 1) * 0.5f;
		float v = 1.0f - ((y + 1) * 0.5f);
		geo.positions.emplace_back(x, 0, y);
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
	geo.texcoords.emplace_back(0, 0);
	geo.texcoords.emplace_back(1, 0);
	geo.texcoords.emplace_back(1, 1);
	geo.texcoords.emplace_back(0, 1);
	geo.normals.assign(4, {0, 0, 1});
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
	geo.texcoords.emplace_back(0, 0);
	geo.texcoords.emplace_back(1, 0);
	geo.texcoords.emplace_back(1, 1);
	geo.texcoords.emplace_back(0, 1);
	geo.normals.assign(4, {0, 1, 0});
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
	geo.texcoords.emplace_back(0, 0);
	geo.texcoords.emplace_back(1, 0);
	geo.texcoords.emplace_back(1, 1);
	geo.texcoords.emplace_back(0, 1);
	geo.normals.assign(4, {1, 0, 0});
	return geo;
}

vec3 extents(const Geometry& geo)
{
	vec3 lower, upper;
	bounds(geo, lower, upper);
	return upper - lower;
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
	}
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

} // namespace gengeo
