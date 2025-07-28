//! This module contains the primitive geometry types.
//!
//! We define primitive geometry as any geometry of O(1) complexity, such as points, lines,
//! triangles, and axis-aligned bounding boxes.
//!
//! Notably, we exclude notions like meshes and paths from this module,
//! as they are of variable complexity.

#![allow(dead_code)]
pub mod aabb;
pub mod infinite_triangle;
pub mod ray;
pub mod segment;
pub mod surface_point;
pub mod triangle;
pub mod vector_utils;
