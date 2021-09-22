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
using System;
using System.Runtime.InteropServices;
using UnityEngine;

namespace Source.P2 {
    public static class TetgenParser {
        [DllImport("TetgenForUnity", EntryPoint = "tetrahedralize_unity")]
        public static extern IntPtr Tetrahedralize(string file);

        [DllImport("TetgenForUnity", EntryPoint = "release_unity_memory")]
        public static extern IntPtr ReleaseMemory(IntPtr array);

        /**
         * Parses the given STL data and returns its tetrahedralized representation.
         * This method requires the plugin Tetgen for Unity.
         */
        public static void Tetrahedralize(string text, Transform transform, bool invertYZAxes,
            Vector3 nodeScale, Vector3 nodeOffset,
            out Solid.Node[] nodes, out Solid.Tetrahedron[] tetrahedra, out Solid.Triangle[] triangles) {
            var ptr = Tetrahedralize(text);

            var points = Marshal.ReadInt32(ptr);
            var tetra = Marshal.ReadInt32(ptr, 4);
            var faces = Marshal.ReadInt32(ptr, 8);

            var offset = 12;

            nodes = new Solid.Node[points];
            tetrahedra = new Solid.Tetrahedron[tetra];
            triangles = new Solid.Triangle[faces];

            Debug.Log("POINTS: " + points);
            Debug.Log("FACES: " + faces);
            Debug.Log("TETRA: " + tetra);

            for (var i = 0; i < points; i++) {
                var x = Int32BitsToSingle(Marshal.ReadInt32(ptr, offset));
                var y = Int32BitsToSingle(Marshal.ReadInt32(ptr, offset + 4));
                var z = Int32BitsToSingle(Marshal.ReadInt32(ptr, offset + 8));
                //Debug.Log($"{x}, {y}, {z}");

                var pos = invertYZAxes ? new Vector3(x, z, y) : new Vector3(x, y, z);
                nodes[i] = new Solid.Node(transform.TransformPoint(pos.Mult(nodeScale) + nodeOffset));

                offset += 12;
            }

            for (var i = 0; i < tetra; i++) {
                var a = Marshal.ReadInt32(ptr, offset);
                var b = Marshal.ReadInt32(ptr, offset + 4);
                var c = Marshal.ReadInt32(ptr, offset + 8);
                var d = Marshal.ReadInt32(ptr, offset + 12);
                //Debug.Log($"{a}, {b}, {c}, {d}");
                tetrahedra[i] = new Solid.Tetrahedron(a, b, c, d);
                offset += 16;
            }

            for (var i = 0; i < faces; i++) {
                var a = Marshal.ReadInt32(ptr, offset);
                var b = Marshal.ReadInt32(ptr, offset + 4);
                var c = Marshal.ReadInt32(ptr, offset + 8);
                triangles[i] = new Solid.Triangle(a, b, c);
                //Debug.Log($"{a}, {b}, {c}");
                offset += 12;
            }

            // The pointer is unmanaged. We need to free that memory!
            ReleaseMemory(ptr);
        }

        /**
         * Reinterprets the int as a float. 
         */
        private static float Int32BitsToSingle(int value) {
            return BitConverter.ToSingle(BitConverter.GetBytes(value), 0);
        }
    }
}