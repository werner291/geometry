use nalgebra::{center, Point3, Vector3};

/// Axis-Aligned Bounding Box (AABB) structure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb {
    /// Minimum point of the AABB.
    pub min: Point3<f64>,
    /// Maximum point of the AABB.
    pub max: Point3<f64>,
}

impl Aabb {
    /// Creates an AABB that encompasses all the given points.
    ///
    /// # Arguments
    ///
    /// * `points` - A slice of points to be enclosed by the AABB.
    ///
    /// # Returns
    ///
    /// An AABB that encloses all the given points.
    pub fn around_points(points: impl IntoIterator<Item = Point3<f64>>) -> Self {
        let mut min = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut max = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);

        for pt in points {
            min.x = min.x.min(pt.x);
            min.y = min.y.min(pt.y);
            min.z = min.z.min(pt.z);

            max.x = max.x.max(pt.x);
            max.y = max.y.max(pt.y);
            max.z = max.z.max(pt.z);
        }

        Self { min, max }
    }

    /// Calculates the center point of the AABB.
    ///
    /// # Returns
    ///
    /// The center point of the AABB.
    pub fn center(&self) -> Point3<f64> {
        center(&self.min, &self.max)
    }

    /// Expands the AABB to cover the given point.
    ///
    /// # Arguments
    ///
    /// * `pt` - The point to be covered by the AABB.
    pub fn expand_to_cover(&mut self, pt: Point3<f64>) {
        self.min.x = self.min.x.min(pt.x);
        self.min.y = self.min.y.min(pt.y);
        self.min.z = self.min.z.min(pt.z);

        self.max.x = self.max.x.max(pt.x);
        self.max.y = self.max.y.max(pt.y);
        self.max.z = self.max.z.max(pt.z);
    }

    /// Creates an AABB that represents a single point.
    ///
    /// # Arguments
    ///
    /// * `pt` - The point to be represented by the AABB.
    ///
    /// # Returns
    ///
    /// An AABB that represents the given point.
    pub fn singleton(pt: Point3<f64>) -> Self {
        Self { min: pt, max: pt }
    }

    /// Compute the half-extents of the AABB.
    pub fn half_extents(&self) -> Vector3<f64> {
        (self.max - self.min) / 2.0
    }

    /// Compute the maximum distance of any point in the AABB from a given point.
    pub fn max_squared_distance(&self, from_point: &Point3<f64>) -> f64 {
        ((self.center() - from_point).abs() + self.half_extents()).norm_squared()
    }

    /// Create a version inflated by a given amount.
    pub fn inflated(&self, amount: f64) -> Self {
        Self {
            min: self.min - Vector3::new(amount, amount, amount),
            max: self.max + Vector3::new(amount, amount, amount),
        }
    }
}
