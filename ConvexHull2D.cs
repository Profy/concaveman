using System;
using System.Collections.Generic;
using System.Numerics;

namespace concaveman
{
    // Monotone chain convex hull algorithm
    public static class ConvexHull2D<T> where T : INumber<T>
    {
        private static T Orient2d(T[] p1, T[] p2, T[] p3)
        {
            return ((p2[1] - p1[1]) * (p3[0] - p2[0])) - ((p2[0] - p1[0]) * (p3[1] - p2[1]));
        }

        public static int[] GetConvexHull(T[][] points)
        {
            int[] result;

            int n = points.Length;
            if (n <= 3)
            {
                result = new int[n];
                for (int i = 0; i < n; ++i)
                {
                    result[i] = i;
                }
                return result;
            }

            // Sort point indices
            int[] sorted = new int[n];
            for (int i = 0; i < n; ++i)
            {
                sorted[i] = i;
            }
            Array.Sort(sorted, (a, b) => { return points[a][0] == points[b][0] ? points[a][1].CompareTo(points[b][1]) : (points[a][0] > points[b][0] ? 1 : -1); });

            // Construct upper and lower hulls
            List<int> lower = [sorted[0], sorted[1]];
            List<int> upper = [sorted[0], sorted[1]];
            for (int i = 2; i < n; ++i)
            {
                int idx = sorted[i];
                T[] p = points[idx];

                // Insert into lower list
                int m = lower.Count;
                while (m > 1 && Orient2d(points[lower[m - 2]], points[lower[m - 1]], p) <= T.Zero)
                {
                    m--;
                    lower.RemoveAt(lower.Count - 1);
                }
                lower.Add(idx);

                // Insert into upper list
                m = upper.Count;
                while (m > 1 && Orient2d(points[upper[m - 2]], points[upper[m - 1]], p) >= T.Zero)
                {
                    m--;
                    upper.RemoveAt(upper.Count - 1);
                }
                upper.Add(idx);
            }

            // Merge lists together
            result = new int[upper.Count + lower.Count - 2];
            int ptr = 0;
            for (int i = 0, nl = lower.Count; i < nl; ++i)
            {
                result[ptr++] = lower[i];
            }
            for (int i = upper.Count - 2; i > 0; --i)
            {
                result[ptr++] = upper[i];
            }

            return result;
        }
    }
}
