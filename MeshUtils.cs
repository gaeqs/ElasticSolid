using UnityEngine;

namespace Source.P2 {
    public static class MeshUtils {
        public static Mesh Copy(this Mesh mesh) {
            var newMesh = new Mesh {
                vertices = mesh.vertices,
                triangles = mesh.triangles,
                uv = mesh.uv,
                normals = mesh.normals,
                colors = mesh.colors,
                tangents = mesh.tangents
            };
            return newMesh;
        }
    }
}