#if UNITY_EDITOR
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

using UnityEngine;

public class FasthullTest : MonoBehaviour
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

    [SerializeField] private int m_count = 10000;

    private List<Vector2> _points = new List<Vector2>();

    public void Generate()
    {
        _points.Clear();
        for (int i = 0; i < m_count; i++)
        {
            _points.Add(new Vector2(UnityEngine.Random.Range(0, 1000), UnityEngine.Random.Range(0, 1000)));
        }
    }

    private void OnDrawGizmos()
    {
        if (_points.Count < 3)
        {
            return;
        }

        float[] points = new float[_points.Count * 2];
        for (int i = 0; i < _points.Count; i++)
        {
            points[i << 1] = _points[i].x;
            points[(i << 1) + 1] = _points[i].y;
        }

        // Convex
        Gizmos.color = Color.green;
        int[] convex_out = new int[points.Length];
        int convex_lenght = ConvexHull2D(points, points.Length / 2, convex_out);
        Array.Resize(ref convex_out, convex_lenght);
        for (int i = 0; i < convex_lenght - 1; i++)
        {
            Gizmos.DrawLine(_points[convex_out[i]], _points[convex_out[i + 1]]);
        }
        Gizmos.DrawLine(_points[convex_out[convex_lenght - 1]], _points[convex_out[0]]);

        // Concave
        Gizmos.color = Color.red;
        float[] concave_out = new float[points.Length];
        int concave_lenght = ConcaveHull2D(points, points.Length / 2, convex_out, convex_lenght, 2, 1, concave_out);
        Array.Resize(ref concave_out, concave_lenght * 2);
        for (int i = 0; i < concave_lenght - 1; i++)
        {
            Gizmos.DrawLine(new Vector2(concave_out[i << 1], concave_out[(i << 1) + 1]), new Vector2(concave_out[(i + 1) << 1], concave_out[((i + 1) << 1) + 1]));
        }
        Gizmos.DrawLine(new Vector2(concave_out[(concave_lenght * 2) - 2], concave_out[(concave_lenght * 2) - 1]), new Vector2(concave_out[0], concave_out[1]));
    }
}
#endif