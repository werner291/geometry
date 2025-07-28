use nalgebra::{Point3, Unit, Vector3};

/// Represents a ray in 3D space, defined by an origin point and a direction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ray {
    pub origin: Point3<f64>,
    pub dir: Unit<Vector3<f64>>,
}

impl Ray {
    /// Creates a new `Ray` from an origin point and a direction vector.
    ///
    /// # Parameters
    /// - `origin`: The origin point of the ray.
    /// - `dir`: The direction vector of the ray.
    ///
    /// # Returns
    /// A `Ray` object representing the ray.
    pub fn new(origin: Point3<f64>, dir: Unit<Vector3<f64>>) -> Self {
        Self { origin, dir }
    }

    /// Returns the point at a given distance along the ray.
    ///
    /// # Parameters
    /// - `t`: The distance along the ray.
    ///
    /// # Returns
    /// The point at the given distance along the ray.
    pub fn at(&self, t: f64) -> Point3<f64> {
        self.origin + self.dir.as_ref() * t
    }

    // #[cfg(test)]
    pub fn geogebra_string(&self) -> String {
        format!(
            "Ray(({:}, {:}, {:}), Vector(({:}, {:}, {:})))",
            self.origin.x, self.origin.y, self.origin.z, self.dir.x, self.dir.y, self.dir.z
        )
    }
}
