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
using UnityEngine;

namespace Source.P2 {
    /**
     * A fixed is responsible of fixing the nodes inside it's bound.
     */
    public class SolidFixer : MonoBehaviour {
        #region UnityVariables

        public UpdatePolicy updatePolicy = UpdatePolicy.Never;
        public Bounds bounds;
        public bool visualizeBounds = true;

        #endregion

        #region PrivateVariables

        private Solid.Node[] _nodes;
        private List<int> _fixedNodes;
        private List<Vector3> _fixedNodesLocalPositions;

        #endregion

        protected void Start() {
            _nodes = GetComponent<Solid>().Nodes;
            _fixedNodes = new List<int>();
            _fixedNodesLocalPositions = new List<Vector3>();
            RefreshFixedStates();
        }

        protected void FixedUpdate() {
            // Refreshes all node's positions. This allows the fixer to be moved along all fixed nodes!o
            for (var i = 0; i < _fixedNodes.Count; i++) {
                _nodes[_fixedNodes[i]].Position = transform.TransformPoint(_fixedNodesLocalPositions[i]);
            }

            if (updatePolicy == UpdatePolicy.Always) {
                RefreshFixedStates();
            }
        }

        /**
         * This method draws a transparent red cube representing the fixer's bounds.
         */
        private void OnDrawGizmosSelected() {
            if (!visualizeBounds) return;
            var center = transform.TransformPoint(bounds.center);
            var size = bounds.size;
            Gizmos.color = new Color(1.0f, 0.0f, 0.0f, 0.5f);
            Gizmos.DrawCube(center, size);
        }

        /**
         * This method searches for nodes to fix.
         * It also unfix any of the previous fixed nodes.
         */
        protected virtual void RefreshFixedStates() {
            _fixedNodes.ForEach(node => _nodes[node].Fixed = false);
            _fixedNodes.Clear();
            _fixedNodesLocalPositions.Clear();

            for (var i = 0; i < _nodes.Length; i++) {
                // Check if the node is already fixed by another fixer.
                if (_nodes[i].Fixed) continue;

                var localPosition = transform.InverseTransformPoint(_nodes[i].Position);
                if (!bounds.Contains(localPosition)) continue;

                // Bounds are in local space.
                _nodes[i].Fixed = true;
                _fixedNodes.Add(i);
                _fixedNodesLocalPositions.Add(localPosition);
            }
        }

        protected void Clear() {
            // If the fixer is disabled, unfix all fixed nodes.
            _fixedNodes.ForEach(index => _nodes[index].Fixed = false);
            _fixedNodes.Clear();
            _fixedNodesLocalPositions.Clear();
        }

        private void OnDisable() {
            Clear();
        }

        private void OnEnable() {
            // Fixes the nodes on enable. OnEnable() is always called before Start() when the
            // game is started. We don't want nodes to be fixed in that time, so we check if the list is null
            // (aka the Start() method was not called yet).
            if (_fixedNodes != null) {
                RefreshFixedStates();
            }
        }

        public enum UpdatePolicy {
            /**
             * Never refresh the fixed nodes.
             */
            Never,

            /**
             * Check the fixed nodes every Update().
             */
            Always,
        }
    }
}