using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace concaveman
{
    internal class Program
    {
        static int[] data = [
            8, 62, 11, 66, 57, 17, 57, 19, 76, 26, 79, 29, 36, 56, 38, 56, 92, 77, 96, 80, 87, 70, 90, 74,
            43, 41, 47, 43, 0, 58, 2, 62, 76, 86, 80, 89, 27, 13, 27, 15, 71, 63, 75, 67, 25, 2, 27, 2, 87,
            6, 88, 6, 22, 90, 23, 93, 22, 89, 22, 93, 57, 11, 61, 13, 61, 55, 63, 56, 17, 85, 21, 87, 33,
            43, 37, 43, 6, 1, 7, 3, 80, 87, 80, 87, 23, 50, 26, 52, 58, 89, 58, 89, 12, 30, 15, 34, 32, 58,
            36, 61, 41, 84, 44, 87, 44, 18, 44, 19, 13, 63, 15, 67, 52, 70, 54, 74, 57, 59, 58, 59, 17, 90,
            20, 92, 48, 53, 52, 56, 92, 68, 92, 72, 26, 52, 30, 52, 56, 23, 57, 26, 88, 48, 88, 48, 66, 13,
            67, 15, 7, 82, 8, 86, 46, 68, 50, 68, 37, 33, 38, 36, 6, 15, 8, 18, 85, 36, 89, 38, 82, 45, 84,
            48, 12, 2, 16, 3, 26, 15, 26, 16, 55, 23, 59, 26, 76, 37, 79, 39, 86, 74, 90, 77, 16, 75, 18,
            78, 44, 18, 45, 21, 52, 67, 54, 71, 59, 78, 62, 78, 24, 5, 24, 8, 64, 80, 64, 83, 66, 55, 70,
            55, 0, 17, 2, 19, 15, 71, 18, 74, 87, 57, 87, 59, 6, 34, 7, 37, 34, 30, 37, 32, 51, 19, 53, 19,
            72, 51, 73, 55, 29, 45, 30, 45, 94, 94, 96, 95, 7, 22, 11, 24, 86, 45, 87, 48, 33, 62, 34, 65,
            18, 10, 21, 14, 64, 66, 67, 67, 64, 25, 65, 28, 27, 4, 31, 6, 84, 4, 85, 5, 48, 80, 50, 81, 1,
            61, 3, 61, 71, 89, 74, 92, 40, 42, 43, 43, 27, 64, 28, 66, 46, 26, 50, 26, 53, 83, 57, 87, 14,
            75, 15, 79, 31, 45, 34, 45, 89, 84, 92, 88, 84, 51, 85, 53, 67, 87, 67, 89, 39, 26, 43, 27, 47,
            61, 47, 63, 23, 49, 25, 53, 12, 3, 14, 5, 16, 50, 19, 53, 63, 80, 64, 84, 22, 63, 22, 64, 26,
            66, 29, 66, 2, 15, 3, 15, 74, 77, 77, 79, 64, 11, 68, 11, 38, 4, 39, 8, 83, 73, 87, 77, 85, 52,
            89, 56, 74, 60, 76, 63, 62, 66, 65, 67
        ];

        static void Main()
        {
            var indexT = new Flatbush<int>(data.Length / 4);

            for (var i = 0; i < data.Length; i += 4)
            {
                indexT.Add(data[i], data[i + 1], data[i + 2], data[i + 3]);
            }
            indexT.Finish();

            var idsT = indexT.Search(40, 40, 60, 60);

            var resultsT = new List<int>();
            for (var i = 0; i < idsT.Count; i++)
            {
                resultsT.Add(data[4 * idsT[i]]);
                resultsT.Add(data[4 * idsT[i] + 1]);
                resultsT.Add(data[4 * idsT[i] + 2]);
                resultsT.Add(data[4 * idsT[i] + 3]);
            }

            //assert.deepEqual(results.sort(compare), [57, 59, 58, 59, 48, 53, 52, 56, 40, 42, 43, 43, 43, 41, 47, 43].sort(compare));



            float[][] points = [[0.8f, 0.2f], [1, 0], [-1, 1], [1, 1], [0.5f, 0.5f], [-0.5f, 0.5f]];




            int[] convex = ConvexHull2D<float>.GetConvexHull(points);
            //float[][] concave = ConcaveHull2D<float>.GetConcaveHull(points, convex, 4, 2, 1);

            //Console.WriteLine(ConvexHull2D<float>.GetConvexHull(points));
            Console.ReadKey();
        }
    }

    /*
    // Monotone chain convex hull algorithm
    public static class ConvexHull2D<T> where T : INumber<T>
    {
        private static T Orient2d(T[] p1, T[] p2, T[] p3)
        {
            return (p2[1] - p1[1]) * (p3[0] - p2[0]) - (p2[0] - p1[0]) * (p3[1] - p2[1]);
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

    // Concaveman concave hull algorithm
    public static class ConcaveHull2D<T> where T : INumber<T>, IMinMaxValue<T>
    {
        private class CompareDist : IComparer<T>
        {
            public int Compare(T x, T y)
            {
                if (x > y)
                {
                    return 1;
                }
                if (x < y)
                {
                    return -1;
                }
                return 0;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="points"></param>
        /// <param name="hull">Start with a convex hull of the points</param>
        /// <param name="maxChildren"></param>
        /// <param name="concavity">A relative measure of concavity; higher value means simpler hull.</param>
        /// <param name="lengthThreshold">When a segment goes below this length threshold, it won't be drilled down further.</param>
        /// <returns></returns>
        public static T[][] GetConcaveHull(T[][] points, int[] hull, int maxChildren, T concavity, T lengthThreshold)
        {
            // exit if hull includes all points already
            if (hull.Length == points.Length)
            {
                return points;
            }

            // index the points with an R-tree
            RTree<T[]> tree = new RTree<T[]>(maxChildren);
            for (int i = 0; i < points.Length; i++)
            {
                tree.Insert(points[i], new T[4] { points[i][0], points[i][1], points[i][0], points[i][1] });
            }

            CircularList circList = new CircularList();
            CircularElement last = null;
            List<CircularElement> queue = new List<CircularElement>();

            // turn the convex hull into a linked list and populate the initial edge queue with the nodes
            for (int i = 0; i < hull.Length; i++)
            {
                tree.Erase(points[i], new T[4] { points[i][0], points[i][1], points[i][0], points[i][1] });
                last = circList.Insert(last, points[i]);
                queue.Add(last);
            }

            // index the segments with an R-tree (for intersection checks)
            RTree<CircularElement> segTree = new RTree<CircularElement>(maxChildren);
            for (int i = 0; i < queue.Count; i++)
            {
                Node node = queue[i].data;
                UpdateBBox(queue[i]);
                segTree.Insert(queue[i], new T[4] { node.minX, node.minY, node.maxX, node.maxY });
            }

            T sqConcavity = concavity * concavity;
            T sqLenThreshold = lengthThreshold * lengthThreshold;

            CircularElement elem;
            // process edges one by one
            while (queue.Count > 0)
            {
                elem = queue[0];
                queue.RemoveAt(0);

                T[] a = elem.prev.data.p;
                T[] b = elem.data.p;
                T[] c = elem.next.data.p;
                T[] d = elem.next.next.data.p;

                // skip the edge if it's already short enough
                T sqLen = GetSqDist(b, c);
                if (sqLen < sqLenThreshold)
                {
                    continue;
                }

                T maxSqLen = sqLen / sqConcavity;

                // find the best connection point for the current edge to flex inward to
                bool ok = false;
                T[] p = FindCandidate(tree, a, b, c, d, maxSqLen, segTree, ref ok);

                // if we found a connection and it satisfies our concavity measure
                if (ok && T.Min(GetSqDist(p, b), GetSqDist(p, c)) <= maxSqLen)
                {
                    // connect the edge endpoints through this point and add 2 new edges to the queue
                    queue.Add(elem);
                    queue.Add(elem.Insert(p));

                    // update point and segment indexes
                    Node node = elem.data;
                    Node next = elem.next.data;

                    tree.Erase(p, [p[0], p[1], p[0], p[1]]);
                    segTree.Erase(elem, [node.minX, node.minY, node.maxX, node.maxY]);

                    UpdateBBox(elem);
                    UpdateBBox(elem.next);

                    segTree.Insert(elem, [node.minX, node.minY, node.maxX, node.maxY]);
                    segTree.Insert(elem.next, [next.minX, next.minY, next.maxX, next.maxY]);
                }
            }

            // convert the resulting hull linked list to an array of points
            List<T[]> concave = new List<T[]>();
            elem = last.next;
            while (true)
            {
                concave.Add(elem.data.p);
                if (elem == last)
                {
                    break;
                }
                elem = elem.next;
            }
            return concave.ToArray();
        }

        private static T[] FindCandidate(RTree<Node> tree, T[] a, T[] b, T[] c, T[] d, T maxDist, RTree<CircularElement> segTree, ref bool ok)
        {
            PriorityQueue<RTree<T[]>, T> queue = new PriorityQueue<RTree<T[]>, T>(new CompareDist());
            Node node = tree.Data;

            // search through the point R-tree with a depth-first search using a priority queue in the order of distance to the edge (b, c)
            ok = false;
            while (true)
            {
                foreach (var child in node)
                {
                    T[] pt = [child.Bounds[0], child.Bounds[1]];
                    T dist = child.IsLeaf ? SqSegDist(pt, b, c) : SqSegBoxDist(b, c, child);
                    if (dist > maxDist)
                    {
                        continue; // skip the node if it's farther than we ever need
                    }

                    queue.Enqueue(tree, dist); // inverser le dist ?
                }

                while (queue.Count > 0 && queue.Peek().IsLeaf)
                {
                    queue.TryDequeue(out RTree<T[]> element, out T priority);

                    // skip all points that are as close to adjacent edges (a,b) and (c,d), and points that would introduce self-intersections when connected
                    T[] pt = [element.Bounds[0], element.Bounds[1]];
                    T d0 = SqSegDist(pt, a, b);
                    T d1 = SqSegDist(pt, c, d);

                    if (priority < d0 && priority < d1 && NoIntersections(b, pt, segTree) && NoIntersections(c, pt, segTree))
                    {
                        ok = true;
                        return element.Data;
                    }
                }

                if (queue.Count == 0)
                {
                    break;
                }

                node = queue.Dequeue();
            }

            return new T[2];
        }


        private static T Orient2d(T[] p1, T[] p2, T[] p3)
        {
            return (p2[1] - p1[1]) * (p3[0] - p2[0]) - (p2[0] - p1[0]) * (p3[1] - p2[1]);
        }

        // check if the edges (p1,q1) and (p2,q2) intersect
        private static bool Intersects(T[] p1, T[] q1, T[] p2, T[] q2)
        {
            return (p1[0] != q2[0] || p1[1] != q2[1]) && (q1[0] != p2[0] || q1[1] != p2[1])
                && (Orient2d(p1, q1, p2) > T.Zero) != (Orient2d(p1, q1, q2) > T.Zero)
                && (Orient2d(p2, q2, p1) > T.Zero) != (Orient2d(p2, q2, q1) > T.Zero);
        }

        // square distance between 2 points
        private static T GetSqDist(T[] p1, T[] p2)
        {
            T dx = p1[0] - p2[0];
            T dy = p1[1] - p2[1];
            return (dx * dx) + (dy * dy);
        }

        // square distance from a point to a segment
        private static T SqSegDist(T[] p, T[] p1, T[] p2)
        {
            T x = p1[0];
            T y = p1[1];
            T dx = p2[0] - x;
            T dy = p2[1] - y;

            if (dx != T.Zero || dy != T.Zero)
            {
                T t = (((p[0] - x) * dx) + ((p[1] - y) * dy)) / ((dx * dx) + (dy * dy));

                if (t > T.One)
                {
                    x = p2[0];
                    y = p2[1];
                }
                else if (t > T.Zero)
                {
                    x += dx * t;
                    y += dy * t;
                }
            }

            dx = p[0] - x;
            dy = p[1] - y;

            return (dx * dx) + (dy * dy);
        }

        // segment to segment distance
        private static T SqSegSegDist(T x0, T y0, T x1, T y1, T x2, T y2, T x3, T y3)
        {
            T ux = x1 - x0;
            T uy = y1 - y0;
            T vx = x3 - x2;
            T vy = y3 - y2;
            T wx = x0 - x2;
            T wy = y0 - y2;
            T a = (ux * ux) + (uy * uy);
            T b = (ux * vx) + (uy * vy);
            T c = (vx * vx) + (vy * vy);
            T d = (ux * wx) + (uy * wy);
            T e = (vx * wx) + (vy * wy);
            T D = (a * c) - (b * b);

            T sc, sN, tc, tN;
            T sD = D;
            T tD = D;

            if (D == T.Zero)
            {
                sN = T.Zero;
                sD = T.One;
                tN = e;
                tD = c;
            }
            else
            {
                sN = (b * e) - (c * d);
                tN = (a * e) - (b * d);
                if (sN < T.Zero)
                {
                    sN = T.Zero;
                    tN = e;
                    tD = c;
                }
                else if (sN > sD)
                {
                    sN = sD;
                    tN = e + b;
                    tD = c;
                }
            }

            if (tN < T.Zero)
            {
                tN = T.Zero;
                if (-d < T.Zero)
                {
                    sN = T.Zero;
                }
                else if (-d > a)
                {
                    sN = sD;
                }
                else
                {
                    sN = -d;
                    sD = a;
                }
            }
            else if (tN > tD)
            {
                tN = tD;
                if (-d + b < T.Zero)
                {
                    sN = T.Zero;
                }
                else if (-d + b > a)
                {
                    sN = sD;
                }
                else
                {
                    sN = -d + b;
                    sD = a;
                }
            }

            sc = sN == T.Zero ? T.Zero : sN / sD;
            tc = tN == T.Zero ? T.Zero : tN / tD;

            T cx = ((T.One - sc) * x0) + (sc * x1);
            T cy = ((T.One - sc) * y0) + (sc * y1);
            T cx2 = ((T.One - tc) * x2) + (tc * x3);
            T cy2 = ((T.One - tc) * y2) + (tc * y3);
            T dx = cx2 - cx;
            T dy = cy2 - cy;

            return (dx * dx) + (dy * dy);
        }

        // square distance from a segment bounding box to the given one
        private static T SqSegBoxDist(T[] a, T[] b, RTree<T[]> bbox)
        {
            if (Inside(a, bbox) || Inside(b, bbox))
            {
                return T.Zero;
            }

            T minX = bbox.Bounds[0];
            T minY = bbox.Bounds[1];
            T maxX = bbox.Bounds[2];
            T maxY = bbox.Bounds[3];

            T d1 = SqSegSegDist(a[0], a[1], b[0], b[1], minX, minY, maxX, minY);
            if (d1 == T.Zero)
            {
                return T.Zero;
            }

            T d2 = SqSegSegDist(a[0], a[1], b[0], b[1], minX, minY, minX, maxY);
            if (d2 == T.Zero)
            {
                return T.Zero;
            }

            T d3 = SqSegSegDist(a[0], a[1], b[0], b[1], maxX, minY, maxX, maxY);
            if (d3 == T.Zero)
            {
                return T.Zero;
            }

            T d4 = SqSegSegDist(a[0], a[1], b[0], b[1], minX, maxY, maxX, maxY);
            if (d4 == T.Zero)
            {
                return T.Zero;
            }

            return T.Min(T.Min(d1, d2), T.Min(d3, d4));
        }

        // check if the point (a,b) is inside the box
        private static bool Inside(T[] a, RTree<T[]> bbox)
        {
            return (a[0] >= bbox.Bounds[0]) && (a[0] <= bbox.Bounds[2]) && (a[1] >= bbox.Bounds[1]) && (a[1] <= bbox.Bounds[3]);
        }

        // check if the edge (a,b) doesn't intersect any other edges
        private static bool NoIntersections(T[] a, T[] b, RTree<CircularElement> segTree)
        {
            T minX = T.Min(a[0], b[0]);
            T minY = T.Min(a[1], b[1]);
            T maxX = T.Max(a[0], b[0]);
            T maxY = T.Max(a[1], b[1]);

            var isect = segTree.Intersection([minX, minY, maxX, maxY]);
            foreach (RTree<CircularElement> ch in isect)
            {
                if (Intersects(ch.Data.data.p, ch.Data.next.data.p, a, b))
                {
                    return false;
                }
            }

            return true;
        }

        // update the bounding box of a node's edge
        private static void UpdateBBox(CircularElement elem)
        {
            elem.data.minX = T.Min(elem.data.p[0], elem.next.data.p[0]);
            elem.data.minY = T.Min(elem.data.p[1], elem.next.data.p[1]);
            elem.data.maxX = T.Max(elem.data.p[0], elem.next.data.p[0]);
            elem.data.maxY = T.Max(elem.data.p[1], elem.next.data.p[1]);
        }

        private class RTree<U>
        {
            private readonly int m_maxChildren;

            private bool m_isLeaf;
            public bool IsLeaf => m_isLeaf;

            private U m_data;
            public U Data => m_data;

            private readonly List<RTree<U>> m_children = new List<RTree<U>>();
            public IEnumerable<RTree<U>> Children => m_children;

            private T[] m_bounds;
            public T[] Bounds => m_bounds;

            public RTree(int maxChildren)
            {
                m_maxChildren = maxChildren;
                m_children = new List<RTree<U>>(maxChildren);

                m_isLeaf = false;
                m_data = default;
                m_bounds = new T[4] { T.MaxValue, T.MaxValue, T.MinValue, T.MinValue };
            }

            private RTree(U data, T[] bounds)
            {
                m_isLeaf = true;
                m_data = data;
                m_bounds = bounds;

                if (bounds[0] > bounds[2] || bounds[1] > bounds[3])
                {
                    throw new ArgumentException("Bounds minima have to be less than maxima");
                }
            }

            public void Insert(U data, T[] bounds)
            {
                if (m_isLeaf)
                {
                    throw new InvalidOperationException("Cannot insert into leaves");
                }

                m_bounds = UpdatedBounds(bounds);
                if (m_children.Count < m_maxChildren)
                {
                    m_children.Add(new RTree<U>(data, bounds));
                    return;
                }

                RTree<U> bestChild = m_children[0];
                T bestVolume = Volume(bestChild.UpdatedBounds(bounds));
                for (int i = 1; i < m_children.Count; i++)
                {
                    T v = Volume(m_children[i].UpdatedBounds(bounds));
                    if (v < bestVolume)
                    {
                        bestVolume = v;
                        bestChild = m_children[i];
                    }
                }
                if (!bestChild.m_isLeaf)
                {
                    bestChild.Insert(data, bounds);
                    return;
                }

                RTree<U> leaf = new RTree<U>(bestChild.Data, bestChild.Bounds);
                bestChild.m_isLeaf = false;
                bestChild.m_data = m_data;
                bestChild.m_children.Add(leaf);
                //bestChild.Insert(data, bounds);
            }

            public void Erase(U data, T[] bounds)
            {
                if (m_isLeaf)
                {
                    throw new InvalidOperationException("Cannot erase from leaves");
                }

                if (!Intersects(bounds))
                {
                    return;
                }

                for (int i = 0; i < m_children.Count; i++)
                {
                    if (!m_children[i].m_isLeaf)
                    {
                        m_children[i].Erase(data, bounds);
                    }
                    else if (m_children[i].m_data.Equals(data) && m_children[i].m_bounds.Equals(bounds))
                    {
                        m_children.RemoveAt(i);
                    }
                }
            }

            private void Intersection(T[] bounds, List<RTree<U>> result)
            {
                if (!Intersects(bounds))
                {
                    return;
                }

                if (m_isLeaf)
                {
                    result.Add(this);
                    return;
                }

                foreach (RTree<U> child in m_children)
                {
                    child.Intersection(bounds, result);
                }
            }

            internal List<RTree<U>> Intersection(T[] bounds)
            {
                List<RTree<U>> results = new List<RTree<U>>();
                Intersection(bounds, results);
                return results;
            }

            public bool Intersects(T[] bounds)
            {
                if (m_bounds[0] > bounds[2] || m_bounds[2] < bounds[0] || m_bounds[1] > bounds[3] || m_bounds[3] < bounds[1])
                {
                    return false;
                }

                return true;
            }

            public T[] UpdatedBounds(T[] childBounds)
            {
                return
                [
                    T.Min(childBounds[0], m_bounds[0]),
                    T.Min(childBounds[1], m_bounds[1]),
                    T.Max(childBounds[2], m_bounds[2]),
                    T.Max(childBounds[3], m_bounds[3]),
                ];
            }

            public static T Volume(T[] bounds)
            {
                T res = T.One;
                res *= bounds[2] - bounds[0];
                res *= bounds[3] - bounds[1];
                return res;
            }
        }

        // Rtree node
        private class Node
        {
            public T[] p;
            public T minX;
            public T minY;
            public T maxX;
            public T maxY;

            public Node()
            {
                p = new T[2];
            }

            public Node(T[] p) : this()
            {
                this.p = p;
            }
        }

        // Rtree element
        private class CircularElement
        {
            public Node data;
            public CircularElement prev;
            public CircularElement next;

            public CircularElement()
            {
                prev = null;
                next = null;
            }

            public CircularElement(T[] p) : this()
            {
                data = new Node(p);
            }

            public CircularElement Insert(T[] p)
            {
                CircularElement elem = new CircularElement(p)
                {
                    prev = this,
                    next = next
                };
                next.prev = elem;
                next = elem;
                return elem;
            }
        }

        private class CircularList
        {
            public CircularElement last;

            public CircularList()
            {
                last = null;
            }

            public CircularElement Insert(CircularElement prev, T[] p)
            {
                CircularElement elem = new CircularElement(p);

                if (prev == null && last != null)
                {
                    throw new InvalidOperationException("Once the list is non-empty you must specify where to insert");
                }

                if (prev == null)
                {
                    elem.prev = elem;
                    elem.next = elem;
                }
                else
                {
                    elem.prev = prev;
                    elem.next = prev.next;
                    prev.next.prev = elem;
                    prev.next = elem;
                }

                last = elem;

                return elem;
            }
        }
    
    }
*/
}
