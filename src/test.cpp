#include "gengeo.hpp"
#include "../FastNoise/FastNoise.h"

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
	
	{
		const int s = 16;
		const int m = 4;
		VoxelGrid voxels(s, s, s);
		box(voxels, {m, m, m}, {s-m, s-m, s-m});
		Geometry voxelized = polygonize(voxels);
		writeObj(voxelized, "voxels.obj");
	}

	{
		FastNoise noise;
		noise.SetNoiseType(FastNoise::Simplex);
		noise.SetFrequency(20);

		const int s = 32;
		VoxelGrid voxels(s, s, s);
		voxels.visit([&](int x, int y, int z) {
			const float n = noise.GetNoise(x, y, z);
			return n > 0 ? 1 : 0;
		});
		Geometry voxelized = polygonize(voxels);
		writeObj(voxelized, "voxelnoise.obj");
	}


	return 0;
}
