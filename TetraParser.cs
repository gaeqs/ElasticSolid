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
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using UnityEngine;

namespace Source.P2 {
    public static class TetraParser {
        private static readonly char[] LineSeparators = {'\n', '\t'};

        /**
         * Generates the solid's nodes using the .node file Tetgen outputs.
         */
        public static Solid.Node[] GenerateNodes(TextAsset asset, Transform transform, bool invertYZAxes,
            Vector3 scale, Vector3 offset) {
            var lines = asset.text.Split(LineSeparators)
                .Where(s => !s.Trim().StartsWith("#")).ToArray();
            var amount = ReadAmount(lines[0]);
            var array = new Solid.Node[amount];

            for (var i = 0; i < array.Length; i++) {
                array[i] = new Solid.Node(
                    transform.TransformPoint(ParseNodePosition(lines[i + 1], invertYZAxes).Mult(scale) + offset)
                );
            }

            return array;
        }

        /**
         * Generates the solid's tetrahedra using the .ele file Tetgen outputs.
         */
        public static Solid.Tetrahedron[] GenerateTetrahedra(TextAsset asset) {
            var lines = asset.text.Split(LineSeparators)
                .Where(s => !s.Trim().StartsWith("#")).ToArray();

            var amount = ReadAmount(lines[0]);
            var array = new Solid.Tetrahedron[amount];

            for (var i = 0; i < array.Length; i++) {
                array[i] = ParseTetrahedron(lines[i + 1]);
            }

            return array;
        }

        /**
         * Generates the solid's triangles using the .face file Tetgen outputs.
         */
        public static Solid.Triangle[] GenerateTriangles(TextAsset asset) {
            var lines = asset.text.Split(LineSeparators)
                .Where(s => !s.Trim().StartsWith("#")).ToArray();

            var amount = ReadAmount(lines[0]);
            var array = new List<Solid.Triangle>();

            for (var i = 0; i < amount; i++) {
                if (ParseTriangle(lines[i + 1], out var triangle)) {
                    array.Add(triangle);
                }
            }

            return array.ToArray();
        }


        private static int ReadAmount(string line) {
            return int.Parse(line.Split(' ')[0], CultureInfo.InvariantCulture);
        }

        private static Vector3 ParseNodePosition(string line, bool invertYZAxes) {
            var data = line.Split(' ').Where(s => s.Length > 0).ToArray();

            return new Vector3(
                float.Parse(data[1], CultureInfo.InvariantCulture),
                float.Parse(data[invertYZAxes ? 3 : 2], CultureInfo.InvariantCulture),
                float.Parse(data[invertYZAxes ? 2 : 3], CultureInfo.InvariantCulture)
            );
        }

        private static Solid.Tetrahedron ParseTetrahedron(string line) {
            var data = line.Split(' ').Where(s => s.Length > 0).ToArray();

            return new Solid.Tetrahedron(
                int.Parse(data[1], CultureInfo.InvariantCulture) - 1,
                int.Parse(data[2], CultureInfo.InvariantCulture) - 1,
                int.Parse(data[3], CultureInfo.InvariantCulture) - 1,
                int.Parse(data[4], CultureInfo.InvariantCulture) - 1
            );
        }

        private static bool ParseTriangle(string line, out Solid.Triangle triangle) {
            var data = line.Split(' ').Where(s => s.Length > 0).ToArray();

            if (int.Parse(data[4], CultureInfo.InvariantCulture) == 0) {
                triangle = new Solid.Triangle();
                return false;
            }

            triangle = new Solid.Triangle(
                int.Parse(data[1], CultureInfo.InvariantCulture) - 1,
                int.Parse(data[2], CultureInfo.InvariantCulture) - 1,
                int.Parse(data[3], CultureInfo.InvariantCulture) - 1
            );

            return true;
        }
    }
}