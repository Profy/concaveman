using System;
using System.Runtime.InteropServices;

/// <summary>
/// Quick creation of a convex and concave shell from a 2d point list.
/// </summary>
public static class FastHull
{
#if UNITY_STANDALONE_WIN
    [DllImport("fasthull64.dll", EntryPoint = "ConvexHull2D", CallingConvention = CallingConvention.Cdecl)]
    private static extern int ConvexHull2D(float[] points, int points_lenght, int[] out_index);
    [DllImport("fasthull64.dll", EntryPoint = "ConcaveHull2D", CallingConvention = CallingConvention.Cdecl)]
    private static extern int ConcaveHull2D(float[] points, int points_lenght, int[] hull_points, int hull_points_lenght, float concavity, float threshold, float[] out_points);
#elif UNITY_STANDALONE_LINUX
    [DllImport("fasthull", EntryPoint = "ConvexHull2D", CallingConvention = CallingConvention.Cdecl)]
    private static extern int ConvexHull2D(float[] points, int points_lenght, int[] out_index);
    [DllImport("fasthull", EntryPoint = "ConcaveHull2D", CallingConvention = CallingConvention.Cdecl)]
    private static extern int ConcaveHull2D(float[] points, int points_lenght, int[] hull_points, int hull_points_lenght, float concavity, float threshold, float[] out_points);
#endif

    /// <summary>
    /// Returns the convex hull. It's not the points, but the indices of the points in the list that are returned.
    /// </summary>
    /// <param name="points">The list of points in the format (x, y, x, y...).</param>
    public static int[] GetConvexHull(float[] points)
    {
        int[] convex_out = new int[points.Length];
        int convex_lenght = ConvexHull2D(points, points.Length / 2, convex_out);
        Array.Resize(ref convex_out, convex_lenght);

        return convex_out;
    }

    /// <summary>
    /// Returns the concave hull. It's the list of the positions of the points (x, y, x, y...) on the hull is returned.
    /// </summary>
    /// <param name="points">The list of points in the format (x, y, x, y...).</param>
    /// <param name="concavity">Is a relative measure of concavity. 1 results in a relatively detailed shape, Infinity results in a convex hull. You can use values lower than 1, but they can produce pretty crazy shapes.</param>
    /// <param name="threshold">When a segment length is under this threshold, it stops being considered for further detalization. Higher values result in simpler shapes.</param>
    public static float[] GetConcaveHull(float[] points, float concavity = 2, float threshold = 0)
    {
        int[] convex_out = new int[points.Length];
        int convex_lenght = ConvexHull2D(points, points.Length / 2, convex_out);
        Array.Resize(ref convex_out, convex_lenght);

        float[] concave_out = new float[points.Length];
        int concave_lenght = ConcaveHull2D(points, points.Length / 2, convex_out, convex_lenght, concavity, threshold, concave_out);
        Array.Resize(ref concave_out, concave_lenght * 2);

        return concave_out;
    }
}