# ElasticSolid

Elastic-Solid script for Unity.

[![YouTube Video](https://img.youtube.com/vi/WEVExe6nK8Y/0.jpg)](https://www.youtube.com/watch?v=WEVExe6nK8Y)

## Introduction

This set of scripts can be used to generate elastic solids in Unity.

The script requires a .node, .ele and .face files that can be generated using TetGen and a
mesh wrapping the mesh you want to use. The script is prepared to do this job automatically using a modified
version of TetGen. Unfortunately, I cannot share the modified TetGen due to licensing issues.

## Features

- Basic forces (Gravity, Spring).
- Volume-based force that prevents the mesh from being deformed.
- Custom TetGen parser.
- Automatic spring generator.
- Substeps.
- Damping force.
- Wing force.
- Penalty collisions (normal and implicit modes).
- Fixers with mobility options.
- Optimized code.
