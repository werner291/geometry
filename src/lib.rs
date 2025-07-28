//! A small collection of geometry utilities and algorithms.
//!
//! Note: there is some overlap with Parry3d here, but we intentionally avoid using Parry3d
//! to avoid introducing a dependency on it in this module.
//!
//! Parry3d may be mentioned in this module, but only for converting to and from Parry3d types.

pub mod intersection;
pub mod mesh;
pub mod mesh_shortest_paths;
pub mod primitive;

#[cfg(test)]
pub mod proptest;

// Re-export commonly used types
pub use primitive::triangle::TriangleEdge;
pub use primitive::triangle::Triangle;
pub use primitive::aabb::Aabb;
pub use primitive::ray::Ray;
pub use primitive::segment::Segment;
pub use primitive::infinite_triangle::InfiniteTriangle;
pub use mesh::trimesh::{TriangleMesh, MeshSurfacePoint, VertexIndex, FaceIndex, IndexTriangle};
pub use mesh::halfedge::{HalfEdgeIndex, HalfEdgeTopology};
pub use mesh::convex_hull::compute_chull; 