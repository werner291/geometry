use crate::primitive::ray::Ray;
use nalgebra::{Point3, Unit, Vector3};

/// A structure representing an infinite triangle defined by an origin point and two direction vectors.
#[derive(Debug, Clone)]
pub struct InfiniteTriangle {
    /// The origin point of the infinite triangle.
    pub origin: Point3<f64>,
    /// The direction vector from the origin to the left side of the triangle.
    pub left_dir: Unit<Vector3<f64>>,
    /// The direction vector from the origin to the right side of the triangle.
    pub right_dir: Unit<Vector3<f64>>,
}

impl InfiniteTriangle {
    /// Creates a new `InfiniteTriangle` with the given origin and direction vectors.
    ///
    /// # Arguments
    ///
    /// * `origin` - The origin point of the infinite triangle.
    /// * `left_dir` - The direction vector from the origin to the left side of the triangle.
    /// * `right_dir` - The direction vector from the origin to the right side of the triangle.
    pub fn new(
        origin: Point3<f64>,
        left_dir: Unit<Vector3<f64>>,
        right_dir: Unit<Vector3<f64>>,
    ) -> Self {
        // Check for NaNs:
        assert!(left_dir.norm_squared().is_finite());
        assert!(right_dir.norm_squared().is_finite());

        // Check that the face isn't degenerate:
        assert!(left_dir.cross(&right_dir).norm_squared() > 1.0e-10);

        Self {
            origin,
            left_dir,
            right_dir,
        }
    }

    /// Computes the normal vector of the infinite triangle; this vector is not normalized.
    ///
    /// # Returns
    ///
    /// A `Vector3<f64>` representing the normal vector of the triangle.
    pub fn normal(&self) -> Vector3<f64> {
        self.left_dir.cross(&self.right_dir)
    }

    /// Creates a `Ray` from the origin of the triangle in the direction of the left side.
    ///
    /// # Returns
    ///
    /// A `Ray` originating from the triangle's origin and pointing in the direction of `left_dir`.
    pub fn left_ray(&self) -> Ray {
        Ray::new(self.origin, self.left_dir)
    }

    /// Creates a `Ray` from the origin of the triangle in the direction of the right side.
    ///
    /// # Returns
    ///
    /// A `Ray` originating from the triangle's origin and pointing in the direction of `right_dir`.
    pub fn right_ray(&self) -> Ray {
        Ray::new(self.origin, self.right_dir)
    }

    /// Checks if a given point lies within the infinite triangle.
    ///
    /// Precondition: The point must lie in the plane of the triangle.
    ///
    /// # Arguments
    ///
    /// * `pt` - A reference to the point to be checked.
    ///
    /// # Returns
    ///
    /// `true` if the point lies within the triangle, `false` otherwise.
    pub fn contains_point(&self, pt: &Point3<f64>) -> bool {
        let delta = pt - self.origin;

        // Check if the point is in the plane of the triangle.
        assert!((self.normal().dot(&delta)).abs() < 1.0e-5);

        let n1 = self.left_dir.cross(&self.normal());
        let n2 = self.normal().cross(&self.right_dir);

        delta.dot(&n1) <= 0.0 && delta.dot(&n2) <= 0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitive::infinite_triangle::InfiniteTriangle;
    use rand::{Rng, SeedableRng};

    #[test]
    fn test_contains_function() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(0);

        for _ in 0..100 {
            let origin = Point3::new(
                rng.random_range(-10.0..10.0),
                rng.random_range(-10.0..10.0),
                rng.random_range(-10.0..10.0),
            );
            let left = Unit::new_normalize(Vector3::new(
                rng.random_range(-1.0..1.0),
                rng.random_range(-1.0..1.0),
                rng.random_range(-1.0..1.0),
            ));
            let right = Unit::new_normalize(Vector3::new(
                rng.random_range(-1.0..1.0),
                rng.random_range(-1.0..1.0),
                rng.random_range(-1.0..1.0),
            ));

            let tri = InfiniteTriangle::new(origin, left, right);

            for _ in 0..10 {
                let u = rng.random_range(0.0..1.0);
                let v = rng.random_range(0.0..1.0);
                let pt = origin + *left * u + *right * v;
                assert!(tri.contains_point(&pt));
            }

            let left_normal = -tri.normal().cross(&left);
            let right_normal = -right.cross(&tri.normal());

            for _ in 0..100 {
                let u = rng.random_range(0.0..1.0);
                let v = rng.random_range(0.0..1.0);
                let pt = origin + *left * u + left_normal * v;
                assert!(!tri.contains_point(&pt));
            }

            for _ in 0..100 {
                let u = rng.random_range(0.0..1.0);
                let v = rng.random_range(0.0..1.0);
                let pt = origin + *right * u + right_normal * v;
                assert!(!tri.contains_point(&pt));
            }
        }
    }
}
