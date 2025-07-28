use crate::mesh::convex_hull::compute_chull;
use crate::mesh::trimesh::{FaceIndex, MeshSurfacePoint, TriangleMesh};
use crate::primitive::triangle::Triangle;
use nalgebra::{Point3, Unit, Vector3};
use ordered_float::OrderedFloat;
use proptest::prelude::*;

/// A Strategy to generate f64 values within the range [-50.0, 50.0].
/// though it is set up to be likely to produce zeros, round numbers, and edge cases.
///
/// We picked 50.0 as that's probably by far the biggest number we'll actually
/// see in our particular project.
pub fn scalar_strategy() -> impl Strategy<Value = f64> {
    prop_oneof![
        prop_oneof![
            Just(0.0),
            Just(1.0),
            Just(-1.0),
            Just(10.0),
            Just(-10.0),
            Just(0.5),
            Just(-0.5),
            Just(std::f64::consts::PI),
            Just(-std::f64::consts::PI),
        ],
        -50.0..=50.0,
    ]
}

pub fn arbitrary_vector() -> impl Strategy<Value = Vector3<f64>> {
    (scalar_strategy(), scalar_strategy(), scalar_strategy())
        .prop_map(|(x, y, z)| Vector3::new(x, y, z))
}

pub fn arbitrary_unit_vector() -> impl Strategy<Value = Unit<Vector3<f64>>> {
    arbitrary_vector()
        .prop_filter("Non-zero vector", |v| v.norm() > 0.0001)
        .prop_map(|v| Unit::new_normalize(v))
}

pub fn arbitrary_point() -> impl Strategy<Value = Point3<f64>> {
    arbitrary_vector().prop_map(|v| Point3::from(v))
}

pub fn arbitrary_triangle() -> impl Strategy<Value = Triangle> {
    (arbitrary_point(), arbitrary_point(), arbitrary_point())
        .prop_map(|(a, b, c)| Triangle::new(a, b, c))
}

pub fn arbitrary_points(max_n: usize) -> impl Strategy<Value = Vec<Point3<f64>>> {
    prop::collection::vec(arbitrary_point(), 1..=max_n)
}

pub fn arbitrary_convex_hull(max_n: usize) -> impl Strategy<Value = TriangleMesh> {
    arbitrary_points(max_n).prop_filter_map("Convex hull", |points| compute_chull(&points).ok())
}

pub fn arbitrary_barycentric() -> impl Strategy<Value = (f64, f64, f64)> {
    prop_oneof![
        // One of the vertices
        Just((1.0, 0.0, 0.0)),
        Just((0.0, 1.0, 0.0)),
        Just((0.0, 0.0, 1.0)),
        // Center of one of the edges
        Just((0.5, 0.5, 0.0)),
        Just((0.0, 0.5, 0.5)),
        Just((0.5, 0.0, 0.5)),
        // Middle:
        Just((1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)),
        // Random on one of the edges
        (0.0..=1.0).prop_map(|x| (x, 1.0 - x, 0.0)),
        (0.0..=1.0).prop_map(|x| (0.0, x, 1.0 - x)),
        (0.0..=1.0).prop_map(|x| (1.0 - x, 0.0, x)),
        // Random
        (0.0..=1.0, 0.0..=1.0).prop_map(|(x, y)| {
            let u = x;
            let y = y * (1.0 - u);
            let z = 1.0 - u - y;
            (u, y, z)
        })
    ]
}

pub fn unit_vector_strategy() -> impl Strategy<Value = Unit<Vector3<f64>>> {
    (
        0.0..std::f64::consts::PI * 2.0,
        -std::f64::consts::FRAC_PI_2..std::f64::consts::FRAC_PI_2,
    )
        .prop_map(|(theta, phi)| {
            let x = theta.sin() * phi.cos();
            let y = theta.sin() * phi.sin();
            let z = theta.cos();
            Unit::new_normalize(Vector3::new(x, y, z))
        })
}

pub fn arbitrary_point_on_chull(chull: TriangleMesh) -> impl Strategy<Value = MeshSurfacePoint> {
    (unit_vector_strategy(), arbitrary_barycentric()).prop_map(move |(dir, bary)| {
        // find the triangle with closest normal:
        let triangle = chull
            .triangles()
            .enumerate()
            .max_by_key(|(_, triangle)| {
                let normal = triangle.unit_normal();
                let dot = normal.dot(&*dir);
                OrderedFloat(dot)
            })
            .expect("Chull should have at least one face.");

        MeshSurfacePoint {
            face_id: FaceIndex(triangle.0),
            point: Point3::from(
                triangle.1.a.coords * bary.0
                    + triangle.1.b.coords * bary.1
                    + triangle.1.c.coords * bary.2,
            ),
        }
    })
}

pub fn arbitrary_point_on_arbitrary_chull(
    max_n: usize,
) -> impl Strategy<Value = (MeshSurfacePoint, TriangleMesh)> {
    arbitrary_convex_hull(max_n).prop_flat_map(|chull| {
        let surface_point = arbitrary_point_on_chull(chull.clone());

        (surface_point, Just(chull))
    })
}

pub fn arbitrary_point_pair_on_arbitrary_chull(
    max_n: usize,
) -> impl Strategy<Value = ([MeshSurfacePoint; 2], TriangleMesh)> {
    arbitrary_convex_hull(max_n).prop_flat_map(|chull| {
        let surface_point = (
            (0..chull.indices.len(), arbitrary_barycentric()),
            (0..chull.indices.len(), arbitrary_barycentric()),
        )
            .prop_map({
                let chull = chull.clone();
                move |pts| {
                    [pts.0, pts.1].map(|(face_id, barycentric)| {
                        let triangle = chull.triangle(FaceIndex(face_id));
                        let a = chull.vertices[triangle.a.0];
                        let b = chull.vertices[triangle.b.0];
                        let c = chull.vertices[triangle.c.0];
                        let point = a.coords * barycentric.0
                            + b.coords * barycentric.1
                            + c.coords * barycentric.2;
                        MeshSurfacePoint {
                            face_id: FaceIndex(face_id as usize),
                            point: Point3::from(point),
                        }
                    })
                }
            });

        (surface_point, Just(chull))
    })
}
