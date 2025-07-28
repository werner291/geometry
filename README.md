# Geometry Crate

A small collection of geometry utilities and algorithms for 3D computational geometry.

## Features

- **Primitive Geometry**: Points, lines, triangles, rays, segments, and axis-aligned bounding boxes
- **Mesh Operations**: Triangle mesh data structures, convex hull computation, half-edge topology
- **Shortest Paths**: Geodesic shortest path algorithms on triangular meshes
- **Intersection Tests**: Ray-segment intersections and parametric line intersections
- **Vector Utilities**: Perpendicular vector computation and angle-based vector operations

## Usage

Add this crate to your `Cargo.toml`:

```toml
[dependencies]
geometry = { git = "https://github.com/your-username/geometry" }
```

## Example

```rust
use geometry::{Triangle, TriangleMesh, Point3};
use nalgebra::Point3;

// Create a triangle
let triangle = Triangle::new(
    Point3::new(0.0, 0.0, 0.0),
    Point3::new(1.0, 0.0, 0.0),
    Point3::new(0.0, 1.0, 0.0)
);

// Compute the area
let area = triangle.area();
println!("Triangle area: {}", area);
```

## License

MIT-0 