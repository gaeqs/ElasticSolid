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