using System;
using System.Collections.Generic;
using System.Numerics;

namespace concaveman
{
    internal class Flatbush<T> where T : INumber<T>, IMinMaxValue<T>, IBitwiseOperators<T, T, T>, IShiftOperators<T, T, T>, IBinaryInteger<T>
    {
        private readonly T hilbertMaxT = T.MultiplicativeIdentity * T.CreateChecked(0xFFFF);

        private readonly int _numItems;
        private readonly int _nodeSize;
        private readonly List<int> _levelBounds;

        private readonly T[] _boxes;
        private readonly int[] _indices;
        private readonly PriorityQueue<int, T> _queue = new PriorityQueue<int, T>();

        private int _pos;

        public T MinX { get; private set; }
        public T MinY { get; private set; }
        public T MaxX { get; private set; }
        public T MaxY { get; private set; }

        /// <summary>
        /// Create a Flatbush index that will hold a given number of items.
        /// </summary>
        /// <param name="numItems">The fixed number of 2d boxes to be included in the index.</param>
        /// <param name="nodeSize">Size of the tree node, adjust to tune for particular use case performance (16 by default).</param>
        public Flatbush(int numItems, int nodeSize = 16)
        {
            if (numItems <= 0)
            {
                throw new ArgumentException("numItems must be greater than zero", nameof(numItems));
            }

            _numItems = numItems;
            _nodeSize = Math.Min(Math.Max(nodeSize, 2), 0xFFFF); // hilbertMax

            // calculate the total number of nodes in the R-tree to allocate space for
            // and the index of each tree level (used in search later)
            int n = numItems;
            int numNodes = n;
            _levelBounds = [n * 4];
            do
            {
                n = (int)Math.Ceiling((double)n / _nodeSize);
                numNodes += n;
                _levelBounds.Add(numNodes * 4);
            } while (n != 1);

            _boxes = new T[numNodes * 4];
            _indices = new int[numNodes];
            _pos = 0;

            MinX = T.MaxValue;
            MinY = T.MaxValue;
            MaxX = T.MinValue;
            MaxY = T.MinValue;
        }

        /// <summary>
        /// Add a given rectangle to the index.
        /// </summary>
        /// <returns>A zero-based, incremental number that represents the newly added rectangle.</returns>
        public int Add(T minX, T minY, T maxX, T maxY)
        {
            int index = _pos >> 2;
            T[] boxes = _boxes;
            _indices[index] = index;
            boxes[_pos++] = minX;
            boxes[_pos++] = minY;
            boxes[_pos++] = maxX;
            boxes[_pos++] = maxY;

            if (minX < MinX)
            {
                MinX = minX;
            }
            if (minY < MinY)
            {
                MinY = minY;
            }
            if (maxX > MaxX)
            {
                MaxX = maxX;
            }
            if (maxY > MaxY)
            {
                MaxY = maxY;
            }

            return index;
        }

        /// <summary>
        /// Method to perform the indexing, to be called after adding all the boxes via <see cref="Add"/>.
        /// </summary>
        public void Finish()
        {
            if (_pos >> 2 != _numItems)
            {
                throw new InvalidOperationException($"Added {_pos >> 2} items when expected {_numItems}.");
            }

            // Only one node, skip sorting and just fill the root box
            if (_numItems <= _nodeSize)
            {
                _boxes[_pos++] = MinX;
                _boxes[_pos++] = MinY;
                _boxes[_pos++] = MaxX;
                _boxes[_pos++] = MaxY;
                return;
            }

            T width = MaxX - MinX;
            T height = MaxY - MinY;
            uint[] hilbertValues = new uint[_numItems];

            // map item centers into Hilbert coordinate space and calculate Hilbert values
            int pos = 0;
            for (int i = 0; i < _numItems; i++)
            {
                T minX = _boxes[pos++];
                T minY = _boxes[pos++];
                T maxX = _boxes[pos++];
                T maxY = _boxes[pos++];

                T x = hilbertMaxT * (((minX + maxX) >> 1) - MinX) / width;
                x -= x % T.One; // Floor
                T y = hilbertMaxT * (((minY + maxY) >> 1) - MinY) / height;
                y -= y % T.One; // Floor
                hilbertValues[i] = Hilbert(uint.CreateChecked(x), uint.CreateChecked(y));
            }

            // sort items by their Hilbert value (for packing later)
            Sort(hilbertValues, _boxes, _indices, 0, _numItems - 1, _nodeSize);

            // generate nodes at each tree level, bottom-up
            pos = 0;
            for (int i = 0; i < _levelBounds.Count - 1; i++)
            {
                // generate a parent node for each block of consecutive <nodeSize> nodes
                int end = _levelBounds[i];
                while (pos < end)
                {
                    int nodeIndex = pos;

                    // calculate bbox for the new node
                    T nodeMinX = _boxes[pos++];
                    T nodeMinY = _boxes[pos++];
                    T nodeMaxX = _boxes[pos++];
                    T nodeMaxY = _boxes[pos++];
                    for (int j = 1; j < _nodeSize && pos < end; j++)
                    {
                        nodeMinX = T.Min(nodeMinX, _boxes[pos++]);
                        nodeMinY = T.Min(nodeMinY, _boxes[pos++]);
                        nodeMaxX = T.Max(nodeMaxX, _boxes[pos++]);
                        nodeMaxY = T.Max(nodeMaxY, _boxes[pos++]);
                    }

                    // add the new node to the tree data
                    _indices[_pos >> 2] = nodeIndex;
                    _boxes[_pos++] = nodeMinX;
                    _boxes[_pos++] = nodeMinY;
                    _boxes[_pos++] = nodeMaxX;
                    _boxes[_pos++] = nodeMaxY;
                }
            }
        }

        /// <summary>
        /// Search the index by a bounding box.
        /// </summary>
        /// <param name=""></param>
        /// <param name=""></param>
        /// <param name=""></param>
        /// <param name=""></param>
        /// <param name=""></param>
        /// <returns></returns>
        /// <exception cref="Error"></exception>
        public List<int> Search(T minX, T minY, T maxX, T maxY, Func<int, bool> filterFn = null)
        {
            if (_pos != _boxes.Length)
            {
                throw new InvalidOperationException("Data not yet indexed, call Finish() before Search().");
            }

            int nodeIndex = _boxes.Length - 4;
            Stack<int> stack = new Stack<int>();
            List<int> results = new List<int>();

            while (true)
            {
                int end = Math.Min(nodeIndex + (_nodeSize * 4), UpperBound(nodeIndex, _levelBounds));
                for (int pos = nodeIndex; pos < end; pos += 4)
                {
                    if (maxX < _boxes[pos])
                    {
                        continue;
                    }
                    if (maxY < _boxes[pos + 1])
                    {
                        continue;
                    }
                    if (minX > _boxes[pos + 2])
                    {
                        continue;
                    }
                    if (minY > _boxes[pos + 3])
                    {
                        continue;
                    }

                    int index = _indices[pos >> 2];
                    if (nodeIndex >= _numItems * 4)
                    {
                        stack.Push(index); // Add node to the search queue.

                    }
                    else if (filterFn == null || filterFn(index))
                    {
                        results.Add(index); // Store leaf item.
                    }
                }

                if (stack.Count > 0)
                {
                    nodeIndex = stack.Pop();
                }
                else
                {
                    break;
                }
            }

            return results;
        }

        /// <summary>
        /// Search items in order of distance from the given point.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="maxDistance">Si 0 ou moins alors est egamle à la valeur max de T</param>
        /// <param name="maxResults"></param>
        /// <param name="filterFn"></param>
        /// <returns></returns>
        /// <exception cref="InvalidOperationException"></exception>
        public List<int> Neighbors(T x, T y, T maxDistance, double maxResults = double.MaxValue, Func<int, bool> filterFn = null)
        {
            if (_pos != _boxes.Length)
            {
                throw new InvalidOperationException("Data not yet indexed, call Finish() before Neighbors().");
            }

            if (maxDistance <= T.Zero)
            {
                maxDistance = T.MaxValue;
            }

            int nodeIndex = _boxes.Length - 4;
            List<int> results = new List<int>();
            T maxDistSquared = maxDistance * maxDistance;

            while (nodeIndex > 0)
            {
                int end = Math.Min(nodeIndex + (_nodeSize * 4), UpperBound(nodeIndex, _levelBounds));
                for (int pos = nodeIndex; pos < end; pos += 4)
                {
                    T dx = AxisDist(x, _boxes[pos], _boxes[pos + 2]);
                    T dy = AxisDist(y, _boxes[pos + 1], _boxes[pos + 3]);
                    T dist = (dx * dx) + (dy * dy);
                    if (dist > maxDistSquared)
                    {
                        continue;
                    }

                    int index = _indices[pos >> 2];
                    if (nodeIndex >= _numItems * 4)
                    {
                        _queue.Enqueue(index << 1, dist); // node (use even id)

                    }
                    else if (filterFn == null || filterFn(index))
                    {
                        _queue.Enqueue((index << 1) + 1, dist); // leaf item (use odd id)
                    }
                }

                // pop items from the queue
                while (_queue.Count > 0)
                {
                    _queue.TryPeek(out int elem, out T dist);
                    if ((elem & 1) == 1)
                    {
                        if (dist > maxDistSquared)
                        {
                            _queue.Clear();
                            return results;
                        }

                        results.Add(_queue.Dequeue() >> 1);
                        if (results.Count == maxResults)
                        {
                            _queue.Clear();
                            return results;
                        }
                    }
                }

                nodeIndex = _queue.Count > 0 ? _queue.Dequeue() >> 1 : 0;
            }

            _queue.Clear();
            return results;
        }

        // 1D distance from a value to a range
        private static T AxisDist(T k, T min, T max)
        {
            return k < min ? min - k : k <= max ? T.Zero : k - max;
        }

        // Binary search for the first value in the array bigger than the given.
        private static int UpperBound(int value, List<int> arr)
        {
            int i = 0;
            int j = arr.Count - 1;
            while (i < j)
            {
                int m = (i + j) >> 1;
                if (arr[m] > value)
                {
                    j = m;
                }
                else
                {
                    i = m + 1;
                }
            }
            return arr[i];
        }

        // custom quicksort that partially sorts bbox data alongside the hilbert values
        private static void Sort(uint[] values, T[] boxes, int[] indices, int left, int right, int nodeSize)
        {
            // check against nodeSize (only need to sort down to nodeSize buckets)
            if (left / nodeSize >= right / nodeSize)
            {
                return;
            }

            uint pivot = values[(left + right) >> 1];
            int i = left - 1;
            int j = right + 1;

            while (true)
            {
                do
                {
                    i++;
                }
                while (values[i] < pivot);
                do
                {
                    j--;
                }
                while (values[j] > pivot);
                if (i >= j)
                {
                    break;
                }

                Swap(values, boxes, indices, i, j);
            }

            Sort(values, boxes, indices, left, j, nodeSize);
            Sort(values, boxes, indices, j + 1, right, nodeSize);
        }

        // swap two values and two corresponding boxes
        private static void Swap(uint[] values, T[] boxes, int[] indices, int i, int j)
        {
            (values[j], values[i]) = (values[i], values[j]);
            int k = 4 * i;
            int m = 4 * j;

            T a = boxes[k];
            T b = boxes[k + 1];
            T c = boxes[k + 2];
            T d = boxes[k + 3];
            boxes[k] = boxes[m];
            boxes[k + 1] = boxes[m + 1];
            boxes[k + 2] = boxes[m + 2];
            boxes[k + 3] = boxes[m + 3];
            boxes[m] = a;
            boxes[m + 1] = b;
            boxes[m + 2] = c;
            boxes[m + 3] = d;

            (indices[j], indices[i]) = (indices[i], indices[j]);
        }

        // Fast Hilbert curve algorithm by http://threadlocalmutex.com/. Ported from C++ https://github.com/rawrunprotected/hilbert_curves (public domain)
        private static uint Hilbert(uint x, uint y)
        {
            uint a = x ^ y;
            uint b = 0xFFFF ^ a;
            uint c = 0xFFFF ^ (x | y);
            uint d = x & (y ^ 0xFFFF);

            uint A = a | (b >> 1);
            uint B = (a >> 1) ^ a;
            uint C = (c >> 1) ^ (b & (d >> 1)) ^ c;
            uint D = (a & (c >> 1)) ^ (d >> 1) ^ d;

            a = A; b = B; c = C; d = D;
            A = (a & (a >> 2)) ^ (b & (b >> 2));
            B = (a & (b >> 2)) ^ (b & ((a ^ b) >> 2));
            C ^= (a & (c >> 2)) ^ (b & (d >> 2));
            D ^= (b & (c >> 2)) ^ ((a ^ b) & (d >> 2));

            a = A; b = B; c = C; d = D;
            A = (a & (a >> 4)) ^ (b & (b >> 4));
            B = (a & (b >> 4)) ^ (b & ((a ^ b) >> 4));
            C ^= (a & (c >> 4)) ^ (b & (d >> 4));
            D ^= (b & (c >> 4)) ^ ((a ^ b) & (d >> 4));

            a = A; b = B; c = C; d = D;
            C ^= (a & (c >> 8)) ^ (b & (d >> 8));
            D ^= (b & (c >> 8)) ^ ((a ^ b) & (d >> 8));

            a = C ^ (C >> 1);
            b = D ^ (D >> 1);

            uint i0 = x ^ y;
            uint i1 = b | (0xFFFF ^ (i0 | a));

            i0 = (i0 | (i0 << 8)) & 0x00FF00FF;
            i0 = (i0 | (i0 << 4)) & 0x0F0F0F0F;
            i0 = (i0 | (i0 << 2)) & 0x33333333;
            i0 = (i0 | (i0 << 1)) & 0x55555555;

            i1 = (i1 | (i1 << 8)) & 0x00FF00FF;
            i1 = (i1 | (i1 << 4)) & 0x0F0F0F0F;
            i1 = (i1 | (i1 << 2)) & 0x33333333;
            i1 = (i1 | (i1 << 1)) & 0x55555555;

            return (i1 << 1) | i0;
        }
    }
}