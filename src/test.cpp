#include "gengeo.hpp"
#include "../third-party/FastNoise/FastNoise.h"

using namespace gengeo;

int main(int argc, char* argv[])
{
	Geometry geo = cube();
	scale(geo, 2);
	for (int i = 1; i < 5; ++i) {
		Geometry temp = sphere(i);
		concat(geo, translate(temp, {2.2f * i, 0, 0}));
	}
	Geometry floor = planey();
	concat(geo, scale(translate(floor, {0.5f, -1, 0}), {12, 1, 2}));
	Geometry backwall = planez();
	concat(geo, scale(translate(backwall, {0.5f, 0, -1.1f}), {12, 3, 1}));
	Geometry sidewall = planex();
	concat(geo, translatex(scale(sidewall, 4), -1));

	center(geo);
	writeObj(geo, "test.obj");

	Geometry circle = circley();
	writeObj(circle, "circle.obj");
	
	Geometry b = box({2, 2, 1});
	flipNormals(b);
	writeObj(b, "flippedcube.obj");
	
	if (true)
	{
		const int s = 12;
		const int m = 2;
		VoxelGrid voxels(s, s, s);
		box(voxels, {m, m, m}, {s-m, s-m, s-m});
		Geometry voxelized = polygonize(voxels);
		writeObj(voxelized, "voxels_cube.obj");
		voxelized = polygonizeSmooth(voxels);
		writeObj(voxelized, "marching_cube.obj");
	}

	if (true)
	{
		FastNoise noise;
		noise.SetNoiseType(FastNoise::Simplex);
		noise.SetFrequency(0.1f);

		const int s = 16;
		VoxelGrid voxels(s, s, s);
		voxels.visit([&](int x, int y, int z) {
			const float n = noise.GetNoise(x, y, z);
			return n > 0 ? 1 : 0;
		});
		Geometry voxelized = polygonize(voxels);
		writeObj(voxelized, "voxelnoise.obj");
		voxelized = polygonizeSmooth(voxels);
		writeObj(voxelized, "marching.obj");
	}

	if (true)
	{
		const int s = 12;
		Geometry polygonized = sdf::polygonize(vec3(-s), vec3(s), vec3(0.25f), [s](vec3 pos) {
			float b = sdf::box(pos, vec3(s * 0.6f, s * 0.2f, s * 0.5f));
			float c = sdf::sphere(pos, s * 0.4f);
			return sdf::opUnion(b, c);
		});
		writeObj(polygonized, "sdf.obj");
	}

	return 0;
}
