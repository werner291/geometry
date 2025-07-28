use nalgebra::Point3;

/// A struct representing a line segment in 3D space.
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Segment {
    /// The starting point of the segment.
    pub a: Point3<f64>,
    /// The ending point of the segment.
    pub b: Point3<f64>,
}

impl Segment {
    /// Creates a new `Segment` with the given starting and ending points.
    ///
    /// # Arguments
    ///
    /// * `a` - The starting point of the segment.
    /// * `b` - The ending point of the segment.
    pub fn new(a: Point3<f64>, b: Point3<f64>) -> Self {
        Self { a, b }
    }

    pub fn closest_point(&self, point: Point3<f64>) -> Point3<f64> {
        // Project origin onto the edge.
        let edge_dir = self.b - self.a;
        let origin_left = point - self.a;
        let t = origin_left.dot(&edge_dir) / edge_dir.norm_squared();
        let t = t.clamp(0.0, 1.0); // Clamp to the edge.

        // Compute the point on the edge.
        self.a + edge_dir * t
    }

    pub fn distance_to_point(&self, point: Point3<f64>) -> f64 {
        let closest = self.closest_point(point);
        (closest - point).norm()
    }

    pub fn length(&self) -> f64 {
        (self.b - self.a).norm()
    }

    pub fn geogebra_string(&self) -> String {
        format!(
            "Segment(({:}, {:}, {:}), ({:}, {:}, {:}))",
            self.a.x, self.a.y, self.a.z, self.b.x, self.b.y, self.b.z
        )
    }
}
