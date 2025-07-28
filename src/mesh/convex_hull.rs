use crate::mesh::trimesh::{
    IndexTriangle, TriangleMesh, VertexIndex, dedup_points_vec, deduplicate_vertices,
};
use crate::primitive::triangle::Triangle;
use itertools::Itertools;
use nalgebra::Point3;
use qhull::Qh;

#[derive(Debug)]
pub enum ConvexHullError {
    NotEnoughPoints,
    CoplanarPoints,
}

/// Computes the convex hull of a set of points.
///
/// # Arguments
///
/// * `points` - A slice of `Point3<f64>` representing the points to compute the convex hull for.
///
/// # Returns
///
/// A `TriangleMesh` representing the computed convex hull.
pub fn compute_chull(points: &[Point3<f64>]) -> Result<TriangleMesh, ConvexHullError> {
    let (_, points) = dedup_points_vec(points);

    if points.len() < 3 {
        return Err(ConvexHullError::NotEnoughPoints);
    }

    if points_are_coplanar(&points) {
        return Err(ConvexHullError::CoplanarPoints);
    }

    // Take the average of the points as the center.
    let center = points
        .iter()
        .fold(Point3::origin(), |acc, pt| acc + pt.coords)
        / points.len() as f64;

    let qh = run_qh(&points);

    let vertices: Vec<Point3<f64>> = extract_triangles(&qh)
        // Check that the normal is outside, flip otherwise:
        .map(|triangle| flip_convex(center, triangle))
        .flat_map(|triangle| [triangle.a, triangle.b, triangle.c].into_iter())
        .collect();

    let indices = indices_assuming_triangle_list(&vertices);

    Ok(deduplicate_vertices(&TriangleMesh::new(vertices, indices)))
}

/// Checks if a set of points are coplanar.
///
/// # Arguments
///
/// * `points` - A slice of `Point3<f64>` representing the points to check for coplanarity.
///
/// # Returns
///
/// `true` if the points are coplanar, `false` otherwise.
pub fn points_are_coplanar(points: &[Point3<f64>]) -> bool {
    let a = points[0];
    let b = points[1];

    // Find a point not collinear with the first two points.
    if let Some(c) = points
        .iter()
        .find(|&p| !Triangle::new(a, b, *p).is_degenerate())
    {
        let normal = Triangle::new(a, b, *c).unit_normal();
        let distance = |p: &Point3<f64>| normal.dot(&(p - a));
        points.iter().all(|p| distance(p).abs() < 1.0e-10)
    } else {
        true // Colinear points are a degenerate case of being coplanar.
    }
}

/// Generates a vector of `IndexTriangle` from a list of vertices assuming the vertices are
/// the vertices of a triangle list, in [a,b,c,a,b,c,...] order.
///
/// # Arguments
///
/// * `vertices` - A vector of `Point3<f64>` representing the vertices of the triangles.
///
/// # Returns
///
/// A vector of `IndexTriangle` where each `IndexTriangle` represents a triangle formed by
/// three consecutive vertices in the input vector.
fn indices_assuming_triangle_list(vertices: &[Point3<f64>]) -> Vec<IndexTriangle> {
    (0..vertices.len() / 3)
        .map(|i| {
            IndexTriangle::new(
                VertexIndex(i * 3),
                VertexIndex(i * 3 + 1),
                VertexIndex(i * 3 + 2),
            )
        })
        .collect()
}

/// Flips the triangle if its normal is pointing towards the center.
///
/// # Arguments
///
/// * `center` - The center point to check the normal direction against.
/// * `triangle` - The triangle to potentially flip.
///
/// # Returns
///
/// A `Triangle` that is flipped if its normal was pointing towards the center,
/// otherwise the original triangle.
fn flip_convex(center: Point3<f64>, triangle: Triangle) -> Triangle {
    if triangle.unit_normal().dot(&(triangle.a - center)) < 0.0 {
        Triangle::new(triangle.c, triangle.b, triangle.a)
    } else {
        triangle
    }
}

/// Extracts triangles from the given Qhull object.
///
/// # Arguments
///
/// * `qh` - A reference to a `Qh` object representing the computed convex hull.
///
/// # Returns
///
/// An iterator over `Triangle` objects representing the triangles of the convex hull.
fn extract_triangles<'a>(qh: &'a Qh) -> impl Iterator<Item = Triangle> + 'a {
    qh.faces().flat_map(move |face| {
        let vertices = face
            .vertices()
            .unwrap()
            .iter()
            .map(|v| Point3::new(v.point()[0], v.point()[1], v.point()[2]))
            .collect::<Vec<_>>();

        // Manually triangulate the polygon if it has more than 3 vertices.
        // See https://github.com/LucaCiucci/qhull-rs/issues/15 as qhull-rs seems
        // to ignore the request to triangulate the resulting output.
        if vertices.len() > 3 {
            triangulate_polygon(&vertices)
        } else {
            vec![Triangle::new(vertices[0], vertices[1], vertices[2])]
        }
    })
}

/// Runs the Qhull algorithm on a set of points to compute the convex hull.
///
/// # Arguments
///
/// * `points` - A slice of `Point3<f64>` representing the points to compute the convex hull for.
///
/// # Returns
///
/// A `Qh` object representing the computed convex hull.
///
/// # Panics
///
/// This function will panic if the Qhull algorithm fails to compute the convex hull.
fn run_qh(points: &[Point3<f64>]) -> Qh {
    assert!(
        points.len() >= 3,
        "At least 3 points are required to compute a convex hull"
    );

    Qh::builder()
        .compute(true)
        .triangulate(true)
        .build_from_iter(points.iter().map(|v| [v.x, v.y, v.z]))
        .expect("Failed to compute convex hull")
}

/// Triangulates a polygon by creating a fan of triangles from a central point.
///
/// Notably, this function does not assume the input vertices are sorted in any particular order.
///
/// # Arguments
/// * `vertices` - A vector of vertices representing the polygon to triangulate
///
/// # Returns
/// A vector of triangles that represent the polygon's triangulation
///
/// # Panics
/// - If the input vertex slice has fewer than 3 vertices
///
/// # Notes
/// - The function assumes the vertices form a simple polygon (no self-intersections)
/// - Triangulation is done by sorting vertices around their center point
fn triangulate_polygon(vertices: &[Point3<f64>]) -> Vec<Triangle> {
    // Validate input: need at least 3 vertices to form a polygon
    assert!(vertices.len() >= 3, "Polygon must have at least 3 vertices");

    // Calculate the centroid (average) of all vertices
    let center = vertices
        .iter()
        .fold(Point3::origin(), |acc, pt| acc + pt.coords)
        / vertices.len() as f64;

    // Calculate the plane normal of the first three vertices
    let normal = (vertices[0] - vertices[1])
        .cross(&(vertices[0] - vertices[2]))
        .normalize();

    // Create a perpendicular vector to help with sorting
    let perpendicular = normal.cross(&(vertices[0] - center));

    // Sort vertices by their angle around the center point
    let in_order = vertices
        .iter()
        .map(|v| {
            // Project vertex onto a 2D plane defined by center and first vertex
            let x = (v - center).dot(&perpendicular);
            let y = (v - center).dot(&(vertices[0] - center));
            let angle = y.atan2(x);
            (v, angle)
        })
        // Collect it to be sorted (sorting is a global operation)
        .collect::<Vec<_>>()
        .into_iter()
        // Sort by angle
        .sorted_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        // Extract vertices
        .map(|(v, _)| *v)
        .collect::<Vec<_>>();

    // Create triangles using a fan triangulation method
    // First vertex is the anchor point for all triangles
    let v1 = in_order[0];
    let mut v2 = in_order[1];

    let mut triangles = Vec::new();
    for i in 2..in_order.len() {
        let v3 = in_order[i];
        triangles.push(Triangle::new(v1, v2, v3));
        v2 = v3;
    }
    triangles
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::geometry::proptest::arbitrary_points;
    use proptest::prelude::*;
    use std::collections::HashSet;

    proptest! {
        #[test]
        fn test_compute_chull(points in arbitrary_points(1000)) {
            match compute_chull(&points) {
                Ok(chull) => {

                    print!("{{");
                    for tri in chull.triangles() {
                        print!("Polygon(({}, {}, {}), ({}, {}, {}), ({}, {}, {})),", tri.a.x, tri.a.y, tri.a.z, tri.b.x, tri.b.y, tri.b.z, tri.c.x, tri.c.y, tri.c.z);
                    }
                    println!("}}");

                    // for all vertices of the convex hull, check that it's one of the given points.
                    for vertex in &chull.vertices {
                        assert!(points.iter().any(|pt| (vertex - pt).norm() < FLOAT_COMPARE));
                    }

                    // Check that the triangles aren't degenerate.
                    for tri in chull.triangles() {
                        assert!(tri.unit_normal().norm() > FLOAT_COMPARE);
                    }

                    const FLOAT_COMPARE: f64 = 1.0e-6;

                    // For all points, check that they are inside the convex hull.
                    for pt in points.iter() {
                        for tri in chull.triangles() {
                            assert!(tri.unit_normal().dot(&(*pt - tri.a)) < FLOAT_COMPARE + 1.0e-6);
                            assert!(tri.unit_normal().dot(&(*pt - tri.b)) < FLOAT_COMPARE + 1.0e-6);
                            assert!(tri.unit_normal().dot(&(*pt - tri.c)) < FLOAT_COMPARE + 1.0e-6);
                        }
                    }

                    // Check manifold-ness by ensuring that every edge has a matching reverse edge.
                    let mut edge_matching = HashSet::new();

                    for tri in chull.indices {
                        for edge in tri.edges() {
                            edge_matching.insert(edge);
                        }
                    }

                    // Then, for every edge, check if its counterpart exists:
                    for edge in &edge_matching {
                        assert!(edge_matching.contains(&[edge[1], edge[0]]));
                    }
                },
                Err(ConvexHullError::NotEnoughPoints) => {
                    assert!(dedup_points_vec(&points).0.len() < 3);
                },
                Err(ConvexHullError::CoplanarPoints) => {
                    assert!(points_are_coplanar(&dedup_points_vec(&points).1));
                }
            };
        }
    }
}
