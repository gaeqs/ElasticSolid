using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex32;
using UnityEngine;
using Random = UnityEngine.Random;

namespace Source.P2 {
    [SelectionBase]
    public class Solid : MonoBehaviour {
        #region Unity Variables

        [Header("Files and initialization")]
        [Tooltip(
            "The STL to use. If this file is not present, the nodes file, tetrahedra file and faces file will be used instead. You must rename the file and use the extension .txt.")]
        public TextAsset stlFile;

        [Tooltip("The nodes file. This file will be ignored if the STL file is set.")]
        public TextAsset nodesFile;

        [Tooltip("The tetrahedra file. This file will be ignored if the STL file is set.")]
        public TextAsset tetrasFile;

        [Tooltip("The faces file. This file will be ignored if the STL file is set.")]
        public TextAsset facesFile;

        [Tooltip("Whether the nodes' Y and Z axes should be swapped.")]
        public bool invertYZAxes;

        [Tooltip("The scale of the nodes. This will be applied during initialization.")]
        public Vector3 nodesScale = Vector3.one;

        [Tooltip("The offset of the nodes. This will be applied during initialization.")]
        public Vector3 nodesOffset = Vector3.zero;

        [Header("Basic simulation")] [Tooltip("The integration method to use.")]
        public Integration integration = Integration.Symplectic;

        [Tooltip("The gravity.")] public Vector3 gravity = new Vector3(0, -9.81f, 0);

        [Tooltip(
            "The mass density of each node. The final mass is the product of this density and the volumen assigned to the node.")]
        public float massDensity = 1;

        [Tooltip("The stiffness of each spring. Stabilizes the length of each spring.")]
        public float springStiffness = 100;

        [Tooltip("The stiffness of each tetrahedron. Stabilizes the volume of each tetrahedron.")]
        public float tetraedronStiffness = 100;

        [Tooltip("The amount of steps done in each FixedUpdate call.")]
        public int substeps = 1;

        [Header("Damping")] [Tooltip("Whether the damping force is enable.")]
        public bool enableDampingForce = true;

        [Tooltip("The damping constant value.")]
        public float damping = 1;

        [Header("Wind")] [Tooltip("Whether the wind force is enabled.")]
        public bool enableWindForce;

        [Tooltip("The static wind velocity.")] public Vector3 staticWindVelocity;

        [Tooltip("The frequency of the oscillating wind velocity.")]
        public float windFrequency;

        [Tooltip("The scale of the oscillating wind velocity.")]
        public Vector3 windFrequencyScale;

        [Tooltip("The scale of the noise wind velocity.")]
        public Vector3 windNoiseScale;

        [Tooltip("The friction of this solid with the wind.")]
        public float windFriction;

        [Header("Collision")] [Tooltip("The integration method to use for the collision system.")]
        public PenaltyIntegration penaltyIntegration = PenaltyIntegration.Normal;

        [Tooltip("The penalty force constant value.")]
        public float penaltyForceConstant = 220;

        [Header("Debug")] [Tooltip("Displays the nodes. (Blue)")]
        public bool showNodes;

        [Tooltip("Displays the tetrahedra. (Red)")]
        public bool showTetrahedra;

        [Tooltip("Displays the springs. (Green)")]
        public bool showSprings;

        [Tooltip("Displays the triangles. (Yellow)")]
        public bool showTriangles;

        [Tooltip("Displays the collisions when they happen.")]
        public bool debugCollisions;

        #endregion

        #region Collision variables

        private Vector3 _minPosition, _maxPosition;
        private BoxCollider _collider;
        private HashSet<Collider> _colliders;
        private SphereCollider _sphere;

        #endregion

        #region Properties

        /**
         * This property returns the registered nodes in this solid.
         * This method is used by external components to get and modify the nodes.
         */
        public Node[] Nodes { get; private set; }

        /**
         * This property returns the registered tetrahedra in this solid.
         * This method is used by external components to get and modify the tetrahedra.
         */
        public Tetrahedron[] Tetrahedra { get; private set; }

        /**
         * This property returns the registered triangles in this solid.
         * This method is used by external components to get and modify the triangles.
         */
        public Triangle[] Triangles { get; private set; }

        /**
         * This property returns the registered springs in this solid.
         * This method is used by external components to get and modify the springs.
         */
        public Spring[] Springs { get; private set; }

        /**
         * This property returns the registered vertices in this solid.
         * This method is used by external components to get and modify the vertices.
         */
        public Dictionary<MeshFilter, Vertex[]> Vertices { get; private set; }

        #endregion

        private void Awake() {
            _colliders = new HashSet<Collider>();
            // Add parent colliders to the set.
            foreach (var other in GetComponentsInParent<Collider>()) {
                _colliders.Add(other);
            }

            _collider = gameObject.AddComponent<BoxCollider>();
            _collider.isTrigger = true;
            _minPosition = Vector3.zero;
            _maxPosition = Vector3.zero;

            _sphere = gameObject.AddComponent<SphereCollider>();
            _sphere.radius = 0.0000001f;
            _sphere.isTrigger = true;

            // If the STL file is found, call Tetgen. Else, load the Tetgen output files.
            if (stlFile == null) {
                Nodes = TetraParser.GenerateNodes(nodesFile, transform, invertYZAxes, nodesScale, nodesOffset);
                Tetrahedra = TetraParser.GenerateTetrahedra(tetrasFile);
                Triangles = TetraParser.GenerateTriangles(facesFile);
            }
            else {
                TetgenParser.Tetrahedralize(stlFile.text, transform, invertYZAxes, nodesScale, nodesOffset,
                    out var nodes, out var tetra, out var triangles);
                Nodes = nodes;
                Tetrahedra = tetra;
                Triangles = triangles;
            }

            for (var i = 0; i < Tetrahedra.Length; i++) {
                Tetrahedra[i].InitializeVolume(Nodes, true);
                Tetrahedra[i].RefreshPlanes(Nodes);
            }

            GenerateSprings();
            GenerateVertices();
        }

        #region Update System

        private void Update() {
            // We have to refresh every mesh vertex each render tick so we can see the solid working!
            // Vertices are in world coordinates, while vertex are in local coordinates. That why we have to
            // apply them the inverse transform.
            foreach (var pair in Vertices) {
                var vertices = new Vector3[pair.Value.Length];
                var tr = pair.Key.transform;
                for (var i = 0; i < pair.Value.Length; i++) {
                    vertices[i] = tr.InverseTransformPoint(pair.Value[i].CalculatePosition(Nodes, Tetrahedra));
                }

                pair.Key.mesh.vertices = vertices;
            }
        }

        private void FixedUpdate() {
            // Call step n times. Here we use the Time.fixedDeltaTime and not a fixed value.
            // This value must be divided by the number of steps to give us a correct simulation
            // in terms of time.
            var timePerStep = Time.fixedDeltaTime / substeps;
            for (var step = 0; step < substeps; step++) {
                Step(timePerStep);
            }
        }

        /**
         * The step function. Called by FixedUpdate() every step.
         */
        private void Step(float deltaTime) {
            // Reset the min and max positions. These positions are used by the collision system.
            // Check the collision algorithm explanation comment below for more info.
            _minPosition = _maxPosition = Nodes.Length > 0 ? Nodes[0].Position : Vector3.zero;

            // First we check the collisions.
            UpdateCollisions();

            // Then we apply the gravity, damping and wind forces.
            for (var i = 0; i < Nodes.Length; i++) {
                Nodes[i].ApplyGravityForces(gravity, massDensity);
                if (enableDampingForce) {
                    Nodes[i].ApplyDampingForce(damping);
                }
            }

            // We also apply the wind force if enabled.
            if (enableWindForce) {
                for (var i = 0; i < Triangles.Length; i++) {
                    Triangles[i].ApplyWindForce(Nodes, staticWindVelocity, windFrequency,
                        windFrequencyScale, windNoiseScale, windFriction);
                }
            }

            // After this we apply the spring forces.
            for (var i = 0; i < Springs.Length; i++) {
                Springs[i].ApplySpringForcesToNodes(Nodes, springStiffness, enableDampingForce, damping);
            }

            // We apply the volume forces.
            for (var i = 0; i < Tetrahedra.Length; i++) {
                Tetrahedra[i].ApplyVolumeForces(Nodes, tetraedronStiffness);
            }

            // We solve the new position and velocity.
            // In this step we also use the foreach loop to calculate the minimum and maximum position of the nodes.
            foreach (var node in Nodes) {
                node.Solve(integration, penaltyIntegration, massDensity, deltaTime);
            }


            for (var i = 0; i < Nodes.Length; i++) {
                Nodes[i].Solve(integration, penaltyIntegration, massDensity, deltaTime);
                var pos = Nodes[i].Position;
                if (penaltyIntegration == PenaltyIntegration.None) continue;
                _minPosition = Vector3.Min(_minPosition, pos);
                _maxPosition = Vector3.Max(_maxPosition, pos);
            }

            // Finally, we update the length of the springs.
            for (var i = 0; i < Springs.Length; i++) {
                Springs[i].UpdateLength(Nodes);
            }

            // And the volume of the tetrahedra
            for (var i = 0; i < Tetrahedra.Length; i++) {
                Tetrahedra[i].InitializeVolume(Nodes, false);
            }


            // We update the bounds of the BoxCollider. This is used to filter potential colliders.
            // Check the collision algorithm explanation below for more info about the collision system.
            UpdateColliderBounds();
        }

        #endregion

        #region Collision System

        /*
         * COLLISION ALGORITHM EXPLANATION.
         *
         * Filter pass:
         *
         * Before computing our collisions, we have to filter which colliders in the scene
         * are candidates for a collision. Computing the algorithm with a collider that is
         * far away is pointless and consumes a lot of CPU.
         *
         * That's why we create a small BoxCollider containing the cloth. To do that we need
         * to compute the minimum and maximum points of the current mesh. Now we can let unity
         * handle this primitive collider, using the events OnTriggerEnter(collider) and OnTriggerExit(collider)
         * to listen its collisions with other colliders.
         * 
         * When a collider enters the area, it's added to the collection _colliders.
         * When a collider exits the area, it's removed from such collection.
         *
         * Computation:
         *
         * We call UpdateCollisions() before every step.
         * This method checks for all collisions in _colliders if any node is penetrating it.
         * This is done using the static method Physics.ComputePenetration().
         * Whether this method returns true, it also returns the direction and the distance of the penetration.
         * We use this information to compute our penalty force.
         *
         * Physics.ComputePenetration() requires two colliders. That's why we create at Start()
         * a very small SphereCollider. This SphereCollider allows us to use the method.
         */

        private void OnTriggerEnter(Collider other) {
            _colliders.Add(other);
        }

        private void OnTriggerExit(Collider other) {
            _colliders.Remove(other);
        }

        /**
         * This method updates the collision bounds of the box collider used to filter potential colliders.
         */
        private void UpdateColliderBounds() {
            if (penaltyIntegration == PenaltyIntegration.None) return;
            // We need to transform the min-max bounds system to a center-size bounds system.
            var center = (_maxPosition + _minPosition) / 2;
            var size = _maxPosition - center;
            _collider.center = transform.InverseTransformPoint(center);
            _collider.size = transform.InverseTransformDirection(size * 2);
        }

        /**
         * This method solves the penalty forces of the present collisions.
         */
        private void UpdateCollisions() {
            if (penaltyIntegration == PenaltyIntegration.None) return;
            foreach (var other in _colliders) {
                for (var i = 0; i < Nodes.Length; i++) {
                    var position = Nodes[i].Position;
                    // Small unity hack. Physics.ComputePenetration is VERY fast. We can take advantage
                    // of this method creating a "dummy" sphere collider with a very small radius.
                    // The method returns whether the collision is present, the depth of the penetration
                    // and the direction of such penetration. We can easily use this information to
                    // apply our desired penalty force!
                    if (!Physics.ComputePenetration(
                        other,
                        other.transform.position,
                        other.transform.rotation,
                        _sphere,
                        position,
                        Quaternion.identity,
                        out var direction,
                        out var distance
                    )) continue;

                    if (integration == Integration.Symplectic && penaltyIntegration == PenaltyIntegration.Implicit) {
                        var normalMatrix = direction.ToDenseMatrix();
                        Nodes[i].ImplicitForce -= penaltyForceConstant * (normalMatrix * normalMatrix.Transpose());
                    }
                    else {
                        Nodes[i].Force -= direction * (distance * penaltyForceConstant);
                    }

                    // We can draw the penetration for debug purposes.
                    if (debugCollisions) {
                        Debug.DrawRay(position, direction * distance);
                    }
                }
            }
        }

        #endregion

        #region Debug

        private void OnDrawGizmos() {
            // Check if the solid was generated.
            if (Nodes == null) return;

            // Debug tetrahedra
            if (showTetrahedra) {
                Gizmos.color = Color.red;
                foreach (var tetrahedron in Tetrahedra) {
                    Gizmos.DrawLine(Nodes[tetrahedron.Node1].Position, Nodes[tetrahedron.Node2].Position);
                    Gizmos.DrawLine(Nodes[tetrahedron.Node2].Position, Nodes[tetrahedron.Node3].Position);
                    Gizmos.DrawLine(Nodes[tetrahedron.Node3].Position, Nodes[tetrahedron.Node4].Position);
                    Gizmos.DrawLine(Nodes[tetrahedron.Node4].Position, Nodes[tetrahedron.Node1].Position);
                }
            }

            // Debug nodes
            if (showNodes) {
                Gizmos.color = Color.blue;
                foreach (var node in Nodes) {
                    Gizmos.DrawSphere(node.Position, 0.02f);
                }
            }

            // Debug springs
            if (showSprings) {
                Gizmos.color = Color.green;
                foreach (var spring in Springs) {
                    Gizmos.DrawLine(Nodes[spring.Node1].Position, Nodes[spring.Node2].Position);
                }
            }

            // Debug triangles
            if (showTriangles) {
                Gizmos.color = Color.yellow;
                foreach (var triangle in Triangles) {
                    Gizmos.DrawLine(Nodes[triangle.Node1].Position, Nodes[triangle.Node2].Position);
                    Gizmos.DrawLine(Nodes[triangle.Node2].Position, Nodes[triangle.Node3].Position);
                    Gizmos.DrawLine(Nodes[triangle.Node3].Position, Nodes[triangle.Node1].Position);
                }
            }
        }

        #endregion

        #region Initialization

        /**
         * This method generates the springs.
         *
         * The algorithm used to generate them is the one described in the description of the cloth practice.
         */
        private void GenerateSprings() {
            // First we create the comparer for the snapshots.
            // First we compare the first value. If they are equal, we compare the second value.
            var comparer = Comparer<SpringSnapshot>.Create((o1, o2) => {
                var value = o1.Node1.CompareTo(o2.Node1);
                if (value == 0) value = o1.Node2.CompareTo(o2.Node2);

                // We want to have repeated snapshots in the sorted set. We will merge them later.
                if (value == 0) value = -1;
                return value;
            });

            // We create the set and add the edges. The first item of the snapshot should be the
            // edge with the lowest index. This is required for the algorithm to work.
            // We also distribute the volumen of the tetrahedra here.
            var tree = new SortedSet<SpringSnapshot>(comparer);
            foreach (var t in Tetrahedra) {
                var volume = t.Volume / 6;
                tree.Add(new SpringSnapshot(t.Node1, t.Node2, volume));
                tree.Add(new SpringSnapshot(t.Node2, t.Node3, volume));
                tree.Add(new SpringSnapshot(t.Node3, t.Node4, volume));
                tree.Add(new SpringSnapshot(t.Node4, t.Node1, volume));
                tree.Add(new SpringSnapshot(t.Node1, t.Node3, volume));
                tree.Add(new SpringSnapshot(t.Node2, t.Node4, volume));
            }

            var list = new List<Spring>();

            SpringSnapshot previous = null;
            foreach (var spring in tree) {
                if (previous == null) {
                    previous = spring;
                }
                else {
                    if (previous.Node1 == spring.Node1 && previous.Node2 == spring.Node2) {
                        // If a duplicate is found, discard the new one and add its volume to the old one.
                        previous.Volume += spring.Volume;
                    }
                    else {
                        list.Add(new Spring(
                            previous.Node1,
                            previous.Node2,
                            Nodes[previous.Node1].Distance(Nodes[previous.Node2]),
                            previous.Volume
                        ));
                        previous = spring;
                    }
                }
            }

            // Adds the last value.
            if (previous != null) {
                list.Add(new Spring(
                    previous.Node1,
                    previous.Node2,
                    Nodes[previous.Node1].Distance(Nodes[previous.Node2]),
                    previous.Volume
                ));
            }


            // Finally, we generate the spring array.
            Springs = list.ToArray();
        }

        /**
         * This method search for mesh filters and binds their vertices to this solid.
         */
        private void GenerateVertices() {
            Vertices = new Dictionary<MeshFilter, Vertex[]>();
            foreach (var filter in GetComponentsInChildren<MeshFilter>()) {
                // We have to copy the mesh to make it writable!
                var mesh = filter.mesh.Copy();
                filter.mesh = mesh;

                Vertices[filter] = mesh.vertices
                    .Select(v => filter.transform.TransformPoint(v))
                    .Select(v => new Vertex(v, Nodes, Tetrahedra)).ToArray();
            }
        }

        #endregion

        #region Structures

        /**
         * This structure represents a Node. A node is defined by a position, a velocity and a force.
         * This structure also contains a boolean representing whether the node should be fixed or not.
         */
        public struct Node {
            public Vector3 Position;
            public Vector3 Velocity;
            public Vector3 Force;
            public Matrix<Complex32> ImplicitForce;
            public bool Fixed;
            public float Volume;

            /**
             * Creates the Node at the given Position.
             * The velocity and the force are zero by default.
             * The node is also not fixed by default.
             * <param name="position">The position of the node.</param>
             */
            public Node(Vector3 position) {
                Position = position;
                Velocity = Vector3.zero;
                Force = Vector3.zero;
                ImplicitForce = new DenseMatrix(3);
                Fixed = false;
                Volume = 1;
            }

            /**
             * Returns the distance between this node and the given one.
             * <param name="other">The other node.</param>
             * <returns>The distance.</returns>
             */
            public float Distance(Node other) {
                return (Position - other.Position).magnitude;
            }

            /**
             * Applies the gravity force to this node.
             * 
             * <param name="gravity">The gravity to apply.</param>
             * <param name="massDensity">The mass density of this node.</param>
             */
            public void ApplyGravityForces(Vector3 gravity, float massDensity) {
                Force += gravity * (massDensity * Volume);
            }

            /**
             * Applies the damping force to this node.
             * 
             * <param name="damping">The damping constant.</param>
             */
            public void ApplyDampingForce(float damping) {
                Force -= damping * Velocity;
            }

            /**
             * This method solves the new position and velocity for the node using the current force and
             * the given mass and delta time. The type of integration to apply is given by the integration
             * parameter.
             * <param name="integration">The integration method.</param>
             * <param name="penaltyIntegration">The penalty integration method.</param>
             * <param name="density">The density of this node.</param>
             * <param name="delta">The delta time.</param>
             */
            public void Solve(Integration integration, PenaltyIntegration penaltyIntegration, float density,
                float delta) {
                if (!Fixed) {
                    switch (integration) {
                        case Integration.Explicit:
                            Position += Velocity * delta;
                            SolveVelocity(penaltyIntegration, density * Volume, delta);
                            break;
                        case Integration.Symplectic:
                            SolveVelocity(penaltyIntegration, density * Volume, delta);
                            Position += Velocity * delta;
                            break;
                        default:
                            throw new ArgumentOutOfRangeException(nameof(integration), integration, null);
                    }
                }

                // Reset the forces just after the integration.
                // This allows external code to add forces to this node.
                Force = Vector3.zero;
                ImplicitForce.Clear();
            }

            /**
             * This method solves the velocity for the node.
             */
            private void SolveVelocity(PenaltyIntegration penaltyIntegration, float mass, float delta) {
                // We only calculate the implicit matrix when the implicit penalty method is enabled.
                // This way we save CPU performance!
                if (penaltyIntegration == PenaltyIntegration.Implicit) {
                    var matrix = (DenseMatrix.CreateIdentity(3)
                                  - delta * delta / mass * ImplicitForce).Inverse();
                    var sVel = (Velocity + Force * delta / mass).ToDenseVector();
                    Velocity = (matrix * sVel).ToUnityVector();
                }
                else {
                    Velocity += Force * delta / mass;
                }
            }
        }

        /**
         * This structure is used as a temporal storage in the spring generation.
         */
        public class SpringSnapshot {
            public int Node1;
            public int Node2;
            public float Volume;

            public SpringSnapshot(int nodeA, int nodeB, float volume) {
                Node1 = Math.Min(nodeA, nodeB);
                Node2 = Math.Max(nodeA, nodeB);
                Volume = volume;
            }
        }

        /**
         * This structure represents a spring. A spring is defined by two references to nodes, a current length
         * and an initial length.
         */
        public struct Spring {
            public int Node1;
            public int Node2;
            public float Length;
            public float InitialLength;
            public float Volume;

            /**
             * Create the spring using two nodes references,
             * the initial length and whether this spring is a flexion spring.
             * 
             * <param name="node1">The reference to the first node as an array index.</param>
             * <param name="node2">The reference to the second node as an array index.</param>
             * <param name="length">The initial length.</param>
             * <param name="volume">The volumen assigned to this spring.</param>
             */
            public Spring(int node1, int node2, float length, float volume) {
                Node1 = node1;
                Node2 = node2;
                Length = length;
                InitialLength = length;
                Volume = volume;
            }

            /**
             * Applies the forces related to the spring to the referenced nodes.
             *
             * <param name="nodes">The array of nodes used to get the nodes to modify.</param>
             * <param name="stiffness">The stiffness of the spring.</param>
             * <param name="dampingEnabled">Whether the damping force should be applied.</param>
             * <param name="damping">The damping constant.</param>
             */
            public void ApplySpringForcesToNodes(Node[] nodes, float stiffness, bool dampingEnabled, float damping) {
                var u = nodes[Node1].Position - nodes[Node2].Position;

                //New force formula
                var force = Volume * stiffness * (Length - InitialLength) / (InitialLength * InitialLength * Length) *
                            u;

                if (dampingEnabled) {
                    // Deformation
                    force += damping * Vector3.Dot(u, nodes[Node1].Velocity - nodes[Node2].Velocity) * u;
                }

                nodes[Node1].Force -= force;
                nodes[Node2].Force += force;
            }

            /**
             * Updates the current length of the node.
             * <param name="nodes">The array of nodes used to get the position of the referenced nodes.</param>
             */
            public void UpdateLength(Node[] nodes) {
                Length = nodes[Node1].Distance(nodes[Node2]);
            }
        }

        /**
         * This structure represents a tetrahedron. A tetrahedron is defined by four sorted references to nodes.
         */
        public struct Tetrahedron {
            public int Node1;
            public int Node2;
            public int Node3;
            public int Node4;

            /**
             * The positions of these planes may not be updated. Use RefreshPlanesAndVolume() to update them.
             */
            public readonly Plane[] Planes;

            /**
             * The initial volume of the tetrahedron.
             */
            public float Volume;

            /**
             * The current volume of the tetrahedron.
             */
            public float CurrentVolume;

            /**
             * Creates te tetrahedron using four sored references to nodes.
             * The planes and the volumes are zero by default.
             */
            public Tetrahedron(int node1, int node2, int node3, int node4) {
                Node1 = node1;
                Node2 = node2;
                Node3 = node3;
                Node4 = node4;
                Planes = new Plane[4];
                Volume = CurrentVolume = 0;
            }

            /**
             * Calculates the current volume of this tetrahedron.
             *
             * If the parameter first is true, this method assigns the calculated value to the initial
             * volume and adds it to the nodes.
             */
            public void InitializeVolume(Node[] nodes, bool first) {
                var pos1 = nodes[Node1].Position;
                var pos2 = nodes[Node2].Position;
                var pos3 = nodes[Node3].Position;
                var pos4 = nodes[Node4].Position;
                CurrentVolume = Vector3.Dot(pos4 - pos1, Vector3.Cross(pos2 - pos1, pos3 - pos1));
                if (!first) return;
                Volume = CurrentVolume;

                // Refresh the mass of all nodes:
                var volumePerNode = Volume / 4;
                nodes[Node1].Volume += volumePerNode;
                nodes[Node2].Volume += volumePerNode;
                nodes[Node3].Volume += volumePerNode;
                nodes[Node4].Volume += volumePerNode;
            }

            /**
             * Refreshes the planes delimiting the area of this tetrahedron.
             */
            public void RefreshPlanes(Node[] nodes) {
                var pos1 = nodes[Node1].Position;
                var pos2 = nodes[Node2].Position;
                var pos3 = nodes[Node3].Position;
                var pos4 = nodes[Node4].Position;
                Planes[0] = new Plane(pos1, pos2, pos4);
                Planes[1] = new Plane(pos2, pos3, pos4);
                Planes[2] = new Plane(pos1, pos3, pos2);
                Planes[3] = new Plane(pos4, pos3, pos1);
            }

            /**
             * Returns whether the given position is inside this tetrahedron.
             * Make sure the planes of the tetrahedron are updated before using this method!
             *
             * <param name="position">The position to check.</param>
             */
            public bool IsInside(Vector3 position) {
                return Planes.All(plane => plane.IsInBack(position));
            }

            /**
             * Generates the weights for the given position at the current state of the tetrahedron.
             * Make sure the volume os the tetrahedron is updated!
             */
            public void GenerateWeights(Node[] nodes, Vector3 position, float[] weights) {
                var pos1 = nodes[Node1].Position;
                var pos2 = nodes[Node2].Position;
                var pos3 = nodes[Node3].Position;
                var pos4 = nodes[Node4].Position;

                weights[0] = GenerateVolume(pos2, pos3, position, pos4) / CurrentVolume;
                weights[1] = GenerateVolume(pos1, pos3, pos4, position) / CurrentVolume;
                weights[2] = GenerateVolume(pos1, pos2, position, pos4) / CurrentVolume;
                weights[3] = GenerateVolume(pos1, pos2, pos3, position) / CurrentVolume;
            }

            /**
             * Applies the volume force to each node of the tetrahedron.
             */
            public void ApplyVolumeForces(Node[] nodes, float stiffness) {
                var pos1 = nodes[Node1].Position;
                var pos2 = nodes[Node2].Position;
                var pos3 = nodes[Node3].Position;
                var pos4 = nodes[Node4].Position;
                var constant = stiffness * (CurrentVolume - Volume);
                nodes[Node1].Force += constant * Vector3.Cross(pos3 - pos2, pos4 - pos2).normalized;
                nodes[Node2].Force += constant * Vector3.Cross(pos4 - pos1, pos3 - pos1).normalized;
                nodes[Node3].Force += constant * Vector3.Cross(pos2 - pos1, pos4 - pos1).normalized;
                nodes[Node4].Force += constant * Vector3.Cross(pos3 - pos1, pos2 - pos1).normalized;
            }

            /**
             * Generates the a tetrahedron volume using the given positions.
             * Watch out! Positions must be in order and the volume may be negative!
             */
            private static float GenerateVolume(Vector3 pos1, Vector3 pos2, Vector3 pos3, Vector3 pos4) {
                return Vector3.Dot(pos4 - pos1, Vector3.Cross(pos2 - pos1, pos3 - pos1));
            }
        }

        /**
         * Represents a mesh vertex bind to this solid.
         */
        public readonly struct Vertex {
            public readonly int TetraIndex;
            public readonly Vector3 InitialPosition;
            public readonly float[] Weights;

            /**
             * Creates the vertex and calculates its weights.
             */
            public Vertex(Vector3 position, Node[] nodes, Tetrahedron[] tetrahedra) {
                TetraIndex = -1;
                InitialPosition = position;
                Weights = new float[4];
                for (var i = 0; i < tetrahedra.Length; i++) {
                    if (!tetrahedra[i].IsInside(position)) continue;
                    TetraIndex = i;
                    tetrahedra[i].GenerateWeights(nodes, position, Weights);
                    break;
                }

                if (TetraIndex != -1) return;
                Debug.DrawLine(position, position + Vector3.up / 100, Color.red, 10);
            }

            /**
             * Calculates the current position for this vertex.
             */
            public Vector3 CalculatePosition(Node[] nodes, Tetrahedron[] tetrahedra) {
                if (TetraIndex == -1) {
                    return InitialPosition;
                }

                var t = tetrahedra[TetraIndex];
                return Weights[0] * nodes[t.Node1].Position
                       + Weights[1] * nodes[t.Node2].Position
                       + Weights[2] * nodes[t.Node3].Position
                       + Weights[3] * nodes[t.Node4].Position;
            }
        }

        /**
         * This structure represents a triangle. A triangle. is defined by three references to nodes.
         * Triangles are used to calculate wind forces.
         */
        public struct Triangle {
            public int Node1, Node2, Node3;

            /**
             * Creates the triangle. Indexes must be in order.
             * <param name="node1">The index of the first node.</param>
             * <param name="node2">The index of the second node.</param>
             * <param name="node3">The index of the third node.</param>
             */
            public Triangle(int node1, int node2, int node3) {
                Node1 = node1;
                Node2 = node2;
                Node3 = node3;
            }

            /**
             * Calculates the normal of the triangle.
             * This normal is not normalized. It's magnitude represents the area of the triangle.
             * <param name="nodes">The nodes array.</param>
             * <returns>The normal not normalized.</returns>
             */
            public Vector3 CalculateNormalNotNormalized(Node[] nodes) {
                var p1 = nodes[Node1].Position;
                var p2 = nodes[Node2].Position;
                var p3 = nodes[Node3].Position;
                var ab = p2 - p1;
                var ac = p3 - p1;

                return Vector3.Cross(ab, ac);
            }

            /**
             * Calculates the center of the triangle.
             * <param name="nodes">The nodes array.</param>
             * <returns>The center.</returns>
             */
            public Vector3 CalculateCenter(Node[] nodes) {
                var p1 = nodes[Node1].Position;
                var p2 = nodes[Node2].Position;
                var p3 = nodes[Node3].Position;
                return (p1 + p2 + p3) / 3.0f;
            }

            /**
             * Calculates the velocity of the triangle.
             * <param name="nodes">The nodes array.</param>
             * <returns>The velocity.</returns>
             */
            public Vector3 CalculateVelocity(Node[] nodes) {
                var v1 = nodes[Node1].Velocity;
                var v2 = nodes[Node2].Velocity;
                var v3 = nodes[Node3].Velocity;
                return (v1 + v2 + v3) / 3.0f;
            }

            /**
             * Applies the wind force to this node.
             * This wind implementation has three factors: the static factor,
             * the oscillating factor and the noise factor.
             * 
             * The result velocity is calculating using the next formula:
             * Vw = static + Vo + Vn
             * 
             * Where:
             * Vo = frequencyScale * (sin(frequency * gameTime) + 1) / 2
             * Vn = noiseScale * (randomVector + (1,1,1) * 0.8)
             * 
             * Where: 
             * (x1, y1, z1) * (x2, y2, z2) = (x1 * x2, y1 * y2, z1 * z2)
             * <param name="nodes">The nodes array.</param>
             * <param name="staticWind">The static factor of the wind force.</param>
             * <param name="windFrequency">The frequency of the oscillating factor of the wind force.</param>
             * <param name="windFrequencyScale">The scale of oscillating factor of the wind force.</param>
             * <param name="windNoiseScale">The scale of the noise factor of the wind force.</param>
             * <param name="friction">The friction of the triangle.</param>
             */
            public void ApplyWindForce(Node[] nodes, Vector3 staticWind, float windFrequency,
                Vector3 windFrequencyScale, Vector3 windNoiseScale, float friction) {
                var noise = windNoiseScale.Mult(Random.insideUnitSphere + Vector3.one * 0.8f);
                var windFrequencyForce =
                    windFrequencyScale * (((float) Math.Sin(windFrequency * Time.time) + 1.0f) / 2.0f);

                var velocity = staticWind + noise + windFrequencyForce;

                var normalNotNormalized = CalculateNormalNotNormalized(nodes);
                var area = normalNotNormalized.magnitude / 2.0f;
                var normal = normalNotNormalized.normalized;
                var triangleVelocity = CalculateVelocity(nodes);

                var force = friction * area *
                    Vector3.Dot(normal, velocity - triangleVelocity) * normal / 3.0f;
                nodes[Node1].Force += force;
                nodes[Node2].Force += force;
                nodes[Node3].Force += force;
            }
        }

        #endregion

        /**
         * This enum contains the integration methods supported by this cloth simulation.
         */
        public enum Integration {
            Explicit = 0,
            Symplectic = 1
        }

        /**
         * This enum contains the integration methods supported for the penalty force.
         * The implicit method only works if the general integration method is Sympletic.
         */
        public enum PenaltyIntegration {
            Normal = 0,
            Implicit = 1,
            None = 2
        }
    }
}