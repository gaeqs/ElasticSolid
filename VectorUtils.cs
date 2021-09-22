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
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex32;
using UnityEngine;

namespace Source.P2 {
    public static class VectorUtils {
        public static Vector3 Mult(this Vector3 a, Vector3 b) => new Vector3(a.x * b.x, a.y * b.y, a.z * b.z);

        public static DenseVector ToDenseVector(this Vector3 vector) =>
            new DenseVector(3) {[0] = vector.x, [1] = vector.y, [2] = vector.z};

        public static DenseMatrix ToDenseMatrix(this Vector3 vector) =>
            new DenseMatrix(3, 1) {[0, 0] = vector.x, [1, 0] = vector.y, [2, 0] = vector.z};

        public static Vector3 ToUnityVector(this Vector<Complex32> vector) =>
            new Vector3(vector[0].Real, vector[1].Real, vector[2].Real);
    }
}