using System;
using System.Collections.Generic;
using System.Numerics;

using RBush;

namespace concaveman
{
    // Monotone chain convex hull algorithm.
    public static class ConvexHull2D<T> where T : INumber<T>
    {
        private static T Orient2d(T[] p1, T[] p2, T[] p3)
        {
            return (p2[1] - p1[1]) * (p3[0] - p2[0]) - (p2[0] - p1[0]) * (p3[1] - p2[1]);
        }

        /// <summary>
        /// Returns the convex hull. It's not the points, but the indices of the points in the list that are returned.
        /// </summary>
        /// <param name="points">The list of points in the format {(x, y), (x, y), ...}.</param>
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

            // Sort point indices.
            int[] sorted = new int[n];
            for (int i = 0; i < n; ++i)
            {
                sorted[i] = i;
            }
            Array.Sort(sorted, (a, b) => { return points[a][0] == points[b][0] ? points[a][1].CompareTo(points[b][1]) : (points[a][0] > points[b][0] ? 1 : -1); });

            // Construct upper and lower hulls.
            List<int> lower = [sorted[0], sorted[1]];
            List<int> upper = [sorted[0], sorted[1]];
            for (int i = 2; i < n; ++i)
            {
                int idx = sorted[i];
                T[] p = points[idx];

                // Insert into lower list.
                int m = lower.Count;
                while (m > 1 && Orient2d(points[lower[m - 2]], points[lower[m - 1]], p) <= T.Zero)
                {
                    m--;
                    lower.RemoveAt(lower.Count - 1);
                }
                lower.Add(idx);

                // Insert into upper list.
                m = upper.Count;
                while (m > 1 && Orient2d(points[upper[m - 2]], points[upper[m - 1]], p) >= T.Zero)
                {
                    m--;
                    upper.RemoveAt(upper.Count - 1);
                }
                upper.Add(idx);
            }

            // Merge lists together.
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

    // Concaveman concave hull algorithm.
    public static class ConcaveHull2D<T> where T : INumber<T>, IMinMaxValue<T>, IBitwiseOperators<T, T, T>, IShiftOperators<T, T, T>, IBinaryInteger<T>
    {
        private class Node : ISpatialData
        {
            private readonly T[] m_p;
            public T[] P => m_p;

            public T MinX;
            public T MinY;
            public T MaxX;
            public T MaxY;

            private readonly Envelope m_envelope;
            public ref readonly Envelope Envelope => ref m_envelope;

            private Node m_prev;
            public Node Prev => m_prev;
            private Node m_next;
            public Node Next => m_next;

            public Node(T[] p)
            {
                m_p = p;
                MinX = m_p[0];
                MinY = m_p[1];
                MaxX = m_p[0];
                MaxY = m_p[1];

                m_envelope = new Envelope(double.CreateChecked(p[0]), double.CreateChecked(p[1]), double.CreateChecked(p[0]), double.CreateChecked(p[1]));
            }

            public Node(T[] p, T minX, T minY, T maxX, T maxY)
            {
                m_p = p;
                MinX = minX;
                MinY = minY;
                MaxX = maxX;
                MaxY = maxY;

                m_envelope = new Envelope(double.CreateChecked(minX), double.CreateChecked(minY), double.CreateChecked(maxX), double.CreateChecked(maxY));
            }

            // Create a new point in a doubly linked list.
            public Node(T[] p, Node prev) : this(p)
            {
                if (prev == null)
                {
                    m_prev = this;
                    m_next = this;
                }
                else
                {
                    m_next = prev.m_next;
                    m_prev = prev;
                    prev.m_next.m_prev = this;
                    prev.m_next = this;
                }
            }
        }

        /// <summary>
        /// Returns the concave hull. It's the list of the positions of the points {(x, y), (x, y), ...} on the hull is returned.
        /// </summary>
        /// <param name="points">The list of points in the format {(x, y), (x, y), ...}.</param>
        /// <param name="hull">Concave hull points list.</param>
        /// <param name="lengthThreshold">When a segment length is under this threshold, it stops being considered for further detalization. Higher values result in simpler shapes.</param>
        /// <param name="concavity">Is a relative measure of concavity. 1 results in a relatively detailed shape, Infinity results in a convex hull. You can use values lower than 1, but they can produce pretty crazy shapes.</param>
        public static T[][] GetConcaveHull2D(T[][] points, int[] hull, T lengthThreshold, T concavity)
        {
            // Index the points with an R-tree
            List<Node> pts = new List<Node>(points.Length);
            RBush<Node> tree = new RBush<Node>(16);
            for (int i = 0; i < points.Length; i += 4)
            {
                pts.Add(new Node(points[i]));
            }
            tree.BulkLoad(pts);

            // Turn the convex hull into a linked list and populate the initial edge queue with the nodes
            List<Node> queue = new List<Node>(hull.Length);
            Node last = null;
            for (int i = 0; i < hull.Length; i++)
            {
                T[] p = points[hull[i]];
                tree.Delete(new Node(p));
                last = new Node(p, last);
                queue.Add(last);
            }

            // Index the segments with an R-tree (for intersection checks)
            RBush<Node> segTree = new RBush<Node>(16);
            Node node;
            for (int i = 0; i < queue.Count; i++)
            {
                node = UpdateBBox(queue[i]);
                segTree.Insert(new Node(node.P, node.MinX, node.MinY, node.MaxX, node.MaxY));
            }

            T sqConcavity = concavity * concavity;
            T sqLenThreshold = lengthThreshold * lengthThreshold;

            // Process edges one by one
            while (queue.Count != 0)
            {
                node = queue[0];
                queue.RemoveAt(0);
                T[] a = node.P;
                T[] b = node.Next.P;

                // Skip the edge if it's already short enough
                T sqLen = GetSqDist(a, b);
                if (sqLen < sqLenThreshold)
                {
                    continue;
                }

                T maxSqLen = sqLen / sqConcavity;

                // Find the best connection point for the current edge to flex inward to
                T[] p = FindCandidate(tree, node.Prev.P, a, b, node.Next.Next.P, maxSqLen, segTree);

                // If we found a connection and it satisfies our concavity measure
                if (p != null && T.Min(GetSqDist(p, a), GetSqDist(p, b)) <= maxSqLen)
                {
                    // Connect the edge endpoints through this point and add 2 new edges to the queue
                    queue.Add(node);
                    queue.Add(new Node(p, node));

                    // Update point and segment indexes
                    tree.Delete(new Node(p));
                    segTree.Delete(node);

                    Node n = UpdateBBox(node);
                    segTree.Insert(new Node(n.P, n.MinX, n.MinY, n.MaxX, n.MaxY));
                    n = UpdateBBox(node.Next);
                    segTree.Insert(new Node(n.P, n.MinX, n.MinY, n.MaxX, n.MaxY));
                }
            }

            // Convert the resulting hull linked list to an array of points
            node = last;
            List<T[]> concave = new List<T[]>();
            do
            {
                concave.Add(node.P);
                node = node.Next;
            } while (node != last);
            concave.Add(node.P);

            return concave.ToArray();
        }

        private static T[] FindCandidate(RBush<Node> tree, T[] a, T[] b, T[] c, T[] d, T maxDist, RBush<Node> segTree)
        {
            PriorityQueue<Node, T> queue = new PriorityQueue<Node, T>(new CompareDist());
            RBush<Node>.Node node = tree.Root;

            // Search through the point R-tree with a depth-first search using a priority queue in the order of distance to the edge (b, c).
            while (true)
            {
                for (int i = 0; i < node.Children.Count; i++)
                {
                    Node child = (Node)node.Children[i];
                    T dist = node.IsLeaf ? SqSegDist(child.P, b, c) : SqSegBoxDist(b, c, child);
                    if (dist > maxDist) // Skip the node if it's farther than we ever need.
                    {
                        continue;
                    }
                    queue.Enqueue(child, dist);
                }

                while (queue.Count > 0 && queue.Peek().Next != null)
                {
                    queue.TryDequeue(out Node element, out T priority);

                    // Skip all points that are as close to adjacent edges (a,b) and (c,d), and points that would introduce self-intersections when connected.
                    var d0 = SqSegDist(element.P, a, b);
                    var d1 = SqSegDist(element.P, c, d);
                    if (priority < d0 && priority < d1 && NoIntersections(b, element.P, segTree) && NoIntersections(c, element.P, segTree))
                    {
                        return element.P;
                    }
                }

                if (queue.Count == 0)
                {
                    break;
                }

                // node = queue.Dequeue();
            }

            return null;
        }

        private class CompareDist : IComparer<T>
        {
            public int Compare(T x, T y)
            {
                return x > y ? 1 : x < y ? -1 : 0;
            }
        }

        private static T Orient2d(T[] p1, T[] p2, T[] p3)
        {
            return ((p2[1] - p1[1]) * (p3[0] - p2[0])) - ((p2[0] - p1[0]) * (p3[1] - p2[1]));
        }

        // Square distance from a segment bounding box to the given one.
        private static T SqSegBoxDist(T[] a, T[] b, Node node)
        {
            if (Inside(a, node) || Inside(b, node))
            {
                return T.Zero;
            }

            T minX = node.MinX;
            T minY = node.MinY;
            T maxX = node.MaxX;
            T maxY = node.MaxY;

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
            return d4 == T.Zero ? T.Zero : T.Min(T.Min(d1, d2), T.Min(d3, d4));
        }

        // Check if the point (a,b) is inside the node bbox.
        private static bool Inside(T[] a, Node node)
        {
            return (a[0] >= node.MinX) && (a[0] <= node.MaxX) && (a[1] >= node.MinY) && (a[1] <= node.MaxY);
        }

        // Check if the edge (a,b) doesn't intersect any other edges.
        private static bool NoIntersections(T[] a, T[] b, RBush<Node> segTree)
        {
            T minX = T.Min(a[0], b[0]);
            T minY = T.Min(a[1], b[1]);
            T maxX = T.Max(a[0], b[0]);
            T maxY = T.Max(a[1], b[1]);

            IReadOnlyList<Node> edges = segTree.Search(new Envelope(double.CreateChecked(minX), double.CreateChecked(minY), double.CreateChecked(maxX), double.CreateChecked(maxY)));
            for (var i = 0; i < edges.Count; i++)
            {
                if (Intersects(edges[i].P, edges[i].Next.P, a, b))
                {
                    return false;
                }
            }
            return true;
        }

        // Check if the edges (p1,q1) and (p2,q2) intersect.
        private static bool Intersects(T[] p1, T[] q1, T[] p2, T[] q2)
        {
            return (p1[0] != q2[0] || p1[1] != q2[1]) && (q1[0] != p2[0] || q1[1] != p2[1])
                && (Orient2d(p1, q1, p2) > T.Zero) != (Orient2d(p1, q1, q2) > T.Zero)
                && (Orient2d(p2, q2, p1) > T.Zero) != (Orient2d(p2, q2, q1) > T.Zero);
        }

        // Update the bounding box of a node's edge.
        private static Node UpdateBBox(Node node)
        {
            node.MinX = T.Min(node.P[0], node.Next.P[0]);
            node.MinY = T.Min(node.P[1], node.Next.P[1]);
            node.MaxX = T.Max(node.P[0], node.Next.P[0]);
            node.MaxY = T.Max(node.P[1], node.Next.P[1]);
            return node;
        }

        // Square distance between 2 points.
        private static T GetSqDist(T[] p1, T[] p2)
        {
            T dx = p1[0] - p2[0];
            T dy = p1[1] - p2[1];
            return (dx * dx) + (dy * dy);
        }

        // Square distance from a point to a segment.
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

        // Segment to segment distance.
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
    }
}