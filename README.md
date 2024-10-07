# fasthull

This repository contains an implementation of [Concaveman algorithm](https://github.com/mapbox/concaveman) for extremely fast generation of convex and concave hulls from a list of 2D points.

C++ and C# Wrappers: A Windows DLL and a Linux library written in C++, along with a C# wrapper for the [concaveman-cpp](https://github.com/sadaszewski/concaveman-cpp). This part is based on the work of [kikitte](https://github.com/kikitte/concaveman-csharp). This can be used with Unity and IL2CPP.

.NET Implementation: An implementation of the [Concaveman algorithm](https://github.com/mapbox/concaveman). Using the RBush package and the new INumber<T> interface did not yield satisfactory performance results. 

An implementation of [Flatbush](https://github.com/mourner/flatbush) was developed to replace RBush, but the static constraints it imposes are too restrictive for my needs