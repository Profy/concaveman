#if 0
g++ -std=c++20 -fPIC -O2 -c -o fasthull.o fasthull.cpp
gcc -std=c++20 -fPIC -O2 -shared fasthull.cpp -o libfasthull.so
exit 0
#endif

#if _MSC_VER // this is defined when compiling with Visual Studio
#define EXPORT_API __declspec(dllexport) // Visual Studio needs annotating exported functions with this
#define CDECL __cdecl // Visual Studio needs annotating exported functions with this
#else
#define EXPORT_API // XCode does not need annotating exported functions, so define is empty
#define CDECL // XCode does not need annotating exported functions, so define is empty
#endif

#include "fasthull.h"

extern "C"
{
	EXPORT_API int CDECL ConvexHull2D(float* points, int points_lenght, int* out_index);
	EXPORT_API int CDECL ConcaveHull2D(float* points, int points_lenght, int* hull_points, int hull_points_lenght, float concavity, float threshold, float* out_points);
}

EXPORT_API int CDECL ConvexHull2D(float* points, int points_lenght, int* out_index)
{
	typedef float T;
	typedef std::array<T, 2> point_type;

	std::vector<point_type> p(points_lenght);
	for (auto i = 0; i < points_lenght; i++)
	{
		p[i] = { points[i << 1], points[(i << 1) + 1] };
	}

	auto convex_points = convexHull<T>(p);

	for (auto i = 0; i < convex_points.size(); i++)
	{
		out_index[i] = convex_points[i];
	}

	return static_cast<int>(convex_points.size());
}

EXPORT_API int CDECL ConcaveHull2D(float* points, int points_lenght, int* hull_points, int hull_points_lenght, float concavity, float threshold, float* out_points)
{
	typedef float T;
	typedef std::array<T, 2> point_type;

	std::vector<point_type> p(points_lenght);
	for (auto i = 0; i < points_lenght; i++)
	{
		p[i] = { points[i << 1], points[(i << 1) + 1] };
	}

	std::vector<int> h(hull_points_lenght);
	for (auto i = 0; i < hull_points_lenght; i++)
	{
		h[i] = hull_points[i];
	}

	auto concave_points = concaveHull<T, 16>(p, h, concavity, threshold);

	for (auto i = 0; i < concave_points.size(); i++)
	{
		out_points[i << 1] = concave_points[i][0];
		out_points[(i << 1) + 1] = concave_points[i][1];
	}

	return static_cast<int>(concave_points.size());
}