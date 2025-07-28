use crate::mesh::halfedge::HalfEdgeIndex;
use crate::primitive::ray::Ray;
use crate::primitive::segment::Segment;
use nalgebra::{Point3, Vector3};

const FLOAT_COMPARE: f64 = 1.0e-10;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgePoint {
    /// The edge of the triangle that the point lies on.
    pub edge: HalfEdgeIndex,
    /// The linear interpolation factor between the two vertices of the edge.
    pub t: f64,
    /// The point on the edge.
    pub point: Point3<f64>,
}

/// Computes the intersection parameters of two parametric lines in 2D space.
///
/// # Arguments
///
/// * `pt1` - A point on the first line.
/// * `dir1` - The direction vector of the first line.
/// * `pt2` - A point on the second line.
/// * `dir2` - The direction vector of the second line.
///
/// # Returns
///
/// An array containing the intersection parameters `t1` and `t2` for the first and second lines respectively.
///
/// If the lines are parallel or coincident, `None` is returned.
///
/// # Note
///
/// This function currently only looks at the x and y components of the input vectors;
/// we need to extend it to 3D to handle edge cases.
pub fn parametric_line_intersection(
    pt1: Point3<f64>,
    dir1: Vector3<f64>,
    pt2: Point3<f64>,
    dir2: Vector3<f64>,
) -> Option<[f64; 2]> {
    let cross_dir1_dir2 = dir1.cross(&dir2);
    let delta = pt2 - pt1;

    let det = cross_dir1_dir2.norm_squared();

    if det.abs() < FLOAT_COMPARE {
        // The lines are parallel or coincident.
        return None;
    }

    let t1 = delta.cross(&dir2).dot(&cross_dir1_dir2) / det;
    let t2 = delta.cross(&dir1).dot(&cross_dir1_dir2) / det;

    Some([t1, t2])
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RaySegmentIntersection {
    /// The ray intersects the segment.
    Intersection([f64; 2]),
    /// The ray is parallel to the segment.
    Parallel,
    /// The ray does not intersect the segment.
    NoIntersection,
}

impl RaySegmentIntersection {
    /// Returns the intersection parameters of the ray and the segment if the intersection occurs within the segment bounds.
    ///
    /// # Returns
    ///
    /// An array containing the intersection parameters `t1` and `t2` if the intersection occurs within the segment bounds, otherwise `None`.
    pub fn intersection(&self) -> Option<[f64; 2]> {
        match self {
            RaySegmentIntersection::Intersection(t) => Some(*t),
            _ => None,
        }
    }
}

/// Computes the intersection parameters of a ray and a segment.
///
/// # Arguments
///
/// * `ray` - The ray to intersect with the segment.
/// * `segment` - The segment to intersect with the ray.
///
/// # Returns
///
/// An `Option` containing an array with the intersection parameters `t1` and `t2` if the intersection occurs within the segment bounds, otherwise `None`.
pub fn ray_segment_intersection(ray: Ray, segment: Segment) -> RaySegmentIntersection {
    // TODO: this is completely inappropriate to put this as a general-purpose test,
    //       since the tolerance would be entirely context-dependent.
    let ray_segment_intersection_tolerance: f64 = 1.0e-10;

    match parametric_line_intersection(ray.origin, *ray.dir, segment.a, segment.b - segment.a) {
        None => RaySegmentIntersection::Parallel,
        Some([t1, t2]) => {
            if (-ray_segment_intersection_tolerance..=1.0 + ray_segment_intersection_tolerance)
                .contains(&t2)
            {
                RaySegmentIntersection::Intersection([t1, t2.clamp(0.0, 1.0)])
            } else {
                RaySegmentIntersection::NoIntersection
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::proptest::arbitrary_point;
    use crate::geometry::proptest::arbitrary_unit_vector;
    use proptest::prelude::*;

    proptest! {

        #[test]
        fn test_parametric_line_intersection(
            intersection_point in arbitrary_point(),
            v1 in arbitrary_unit_vector(),
            v2 in arbitrary_unit_vector(),
            min1 in -10.0 .. 0.0,
            max1 in 0.0 .. 10.0,
            min2 in -10.0 .. 0.0,
            max2 in 0.0 .. 10.0,
        ) {
            let segment1 = Segment::new(intersection_point + *v1 * min1, intersection_point + *v1 * max1);
            let segment2 = Segment::new(intersection_point + *v2 * min2, intersection_point + *v2 * max2);

            match parametric_line_intersection(segment1.a, segment1.b - segment1.a, segment2.a, segment2.b - segment2.a) {
                Some([t1, t2]) => {
                    let pt1 = segment1.a + (segment1.b - segment1.a) * t1;
                    let pt2 = segment2.a + (segment2.b - segment2.a) * t2;

                    // Both must be at the intersection point:
                    prop_assert!((pt1 - intersection_point).norm() < 1.0e-6);
                    prop_assert!((pt2 - intersection_point).norm() < 1.0e-6);
                }
                None => {
                    // Check that v1, v2 are parallel or coincident:
                    prop_assert!(
                        (segment1.a - segment2.a).cross(&(segment1.b - segment2.a)).norm() < 1.0e-6
                    );
                }
            }
        }
    }
}
