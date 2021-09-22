/*
MIT License

Copyright (c) 2021 gaeqs

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
using UnityEngine;

namespace Source.P2 {
    /// <summary>
    /// Represents the plane equation p -> ax + by + cz + d = 0;
    ///
    /// This structure was made by me two years ago when I was making a Minecraft clone.
    /// You can check the original code here:
    /// https://github.com/gaeqs/CMine/blob/master/CMineNew/Src/Geometry/Plane.cs
    /// 
    /// </summary>
    public struct Plane {
        private const float DistanceTolerance = 0.00000001f;
        private const float RotationTolerance = 0.00000001f;

        private float _a, _b, _c, _d;


        /// <summary>
        /// Creates a plane using the given plane parameters.
        /// </summary>
        /// <param name="a">The "a" parameter.</param>
        /// <param name="b">The "b" parameter.</param>
        /// <param name="c">The "c" parameter.</param>
        /// <param name="d">The "d" parameter.</param>
        public Plane(float a, float b, float c, float d) {
            _a = a;
            _b = b;
            _c = c;
            _d = d;
        }

        /// <summary>
        /// Creates a plane using a given normal vector and a point in the plane.
        /// </summary>
        /// <param name="normal">The normal vector.</param>
        /// <param name="point">The point in the plane.</param>
        public Plane(Vector3 normal, Vector3 point) {
            normal.Normalize();
            _a = normal.x;
            _b = normal.y;
            _c = normal.z;
            _d = -_a * point.x - _b * point.y - _c * point.z;
        }


        /// <summary>
        /// Creates a plane using the three given points.
        /// All three points will be inside the plane.
        /// </summary>
        /// <param name="point0">The first point.</param>
        /// <param name="point1">The second point.</param>
        /// <param name="point2">The third point.</param>
        public Plane(Vector3 point0, Vector3 point1, Vector3 point2) {
            var normal = Vector3.Cross(point1 - point0, point2 - point0).normalized;
            _a = normal.x;
            _b = normal.y;
            _c = normal.z;
            _d = -_a * point0.x - _b * point0.y - _c * point0.z;
        }

        /// <summary>
        /// The "a" parameter of the equation.
        /// </summary>
        public float A {
            get => _a;
            set => _a = value;
        }

        /// <summary>
        /// The "b" parameter of the equation.
        /// </summary>
        public float B {
            get => _b;
            set => _b = value;
        }


        /// <summary>
        /// The "c" parameter of the equation.
        /// </summary>
        public float C {
            get => _c;
            set => _c = value;
        }


        /// <summary>
        /// THe "d" parameter of the equation.
        /// </summary>
        public float D {
            get => _d;
            set => _d = value;
        }

        /// <summary>
        /// The normal of the plane. This vectors contains the "a", "b" and "c" parameters.
        /// </summary>
        public Vector3 Normal => new Vector3(_a, _b, _c);

        /// <summary>
        /// Normalizes the plane. This is required for a correct distance calculation.
        /// </summary>
        public void Normalize() {
            var length = Mathf.Sqrt(_a * _a + _b * _b + _c * _c);
            _a /= length;
            _b /= length;
            _c /= length;
            _d /= length;
        }

        /**
         * Converts the front of this plane to the back, and the back to the front.
         */
        public void Negate() {
            _a = -_a;
            _b = -_b;
            _c = -_c;
            _d = -_d;
        }


        /// <summary>
        /// Calculates the minimum distance between the plane and the given point.
        /// This distance may be negative. A negative value indicates that the
        /// point is in the back of the plane.
        /// </summary>
        /// <param name="point">The point.</param>
        public float Distance(Vector3 point) {
            return _a * point.x + _b * point.y + _c * point.z + _d;
        }

        /// <summary>
        /// Calculates the minimum distance between the plane and the given point.
        /// This distance may be negative. A negative value indicates that the
        /// point is in the back of the plane.
        /// </summary>
        /// <param name="x">The X value of the point.</param>
        /// <param name="y">The Y value of the point.</param>
        /// <param name="z">The Z value of the point.</param>
        public float Distance(float x, float y, float z) {
            return _a * x + _b * y + _c * z + _d;
        }


        /// <summary>
        /// Returns the minimum distance between the plane and the given point.
        /// </summary>
        /// <param name="point">The point.</param>
        /// <returns></returns>
        public float DistanceAbs(Vector3 point) {
            return Mathf.Abs(Distance(point));
        }

        /// <summary>
        /// Returns the minimum distance between the plane and the given point.
        /// </summary>
        /// <param name="x">The X value of the point.</param>
        /// <param name="y">The Y value of the point.</param>
        /// <param name="z">The Z value of the point.</param>
        /// <returns></returns>
        public float DistanceAbs(float x, float y, float z) {
            return Mathf.Abs(Distance(x, y, z));
        }


        /// <summary>
        /// Returns whether the given point is in the front of the plane.
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        public bool IsInFront(Vector3 point) {
            return Distance(point) > DistanceTolerance;
        }


        /// <summary>
        /// Returns whether the given point is in the front of the plane.
        /// </summary>
        /// <param name="x">The X value of the point.</param>
        /// <param name="y">The Y value of the point.</param>
        /// <param name="z">The Z value of the point.</param>
        /// <returns></returns>
        public bool IsInFront(float x, float y, float z) {
            return Distance(x, y, z) > DistanceTolerance;
        }

        /// <summary>
        /// Returns whether the given point is in the back of the plane.
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        public bool IsInBack(Vector3 point) {
            return Distance(point) < -DistanceTolerance;
        }

        /// <summary>
        /// Returns whether the given point is in the back of the plane.
        /// </summary>
        /// <param name="x">The X value of the point.</param>
        /// <param name="y">The Y value of the point.</param>
        /// <param name="z">The Z value of the point.</param>
        /// <returns></returns>
        public bool IsInBack(float x, float y, float z) {
            return Distance(x, y, z) < -DistanceTolerance;
        }

        /// <summary>
        /// Returns whether the given point is inside the plane.
        /// </summary>
        /// <param name="point">The plane</param>
        /// <returns></returns>
        public bool Intersects(Vector3 point) {
            return DistanceAbs(point) < DistanceTolerance;
        }

        /// <summary>
        /// Returns whether the given plane intersects with this one.
        /// </summary>
        /// <param name="plane">The plane.</param>
        /// <returns></returns>
        public bool Intersects(Plane plane) {
            return Mathf.Abs(plane._d - _d) < DistanceTolerance ||
                   Vector3.Dot(plane.Normal, Normal) < 1 - RotationTolerance;
        }

        /// <summary>
        /// Returns whether the given sphere intersects with the plane.
        /// </summary>
        /// <param name="center">The center of the sphere.</param>
        /// <param name="radius">The radius of the sphere.</param>
        /// <returns></returns>
        public bool IntersectsSphere(Vector3 center, float radius) {
            return DistanceAbs(center) <= radius;
        }
    }
}