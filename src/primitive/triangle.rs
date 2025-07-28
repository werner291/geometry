use crate::primitive::segment::Segment;
use nalgebra::{Point3, Unit, Vector3};
use proptest_derive::Arbitrary;

/// Enum to designate one of the triangle edges.
#[derive(Debug, Clone, Copy, PartialEq, Arbitrary)]
pub enum TriangleEdge {
    EdgeAB,
    EdgeBC,
    EdgeCA,
}

/// A triangle in 3D space defined by three points.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Triangle {
    /// The first vertex of the triangle.
    pub a: Point3<f64>,
    /// The second vertex of the triangle.
    pub b: Point3<f64>,
    /// The third vertex of the triangle.
    pub c: Point3<f64>,
}

impl Triangle {
    pub fn is_degenerate(&self) -> bool {
        // Check if the triangle is degenerate by checking if the vertices are collinear.
        let ab = self.b - self.a;
        let ac = self.c - self.a;
        ab.cross(&ac).norm() < 1e-10
    }
}

impl Triangle {
    pub fn new(a: Point3<f64>, b: Point3<f64>, c: Point3<f64>) -> Self {
        Self { a, b, c }
    }

    pub fn raw_normal(&self) -> Vector3<f64> {
        (self.b - self.a).cross(&(self.c - self.a))
    }

    pub fn unit_normal(&self) -> Unit<Vector3<f64>> {
        Unit::new_normalize(self.raw_normal())
    }

    pub fn edge(&self, edge: TriangleEdge) -> Segment {
        match edge {
            TriangleEdge::EdgeAB => Segment::new(self.a, self.b),
            TriangleEdge::EdgeBC => Segment::new(self.b, self.c),
            TriangleEdge::EdgeCA => Segment::new(self.c, self.a),
        }
    }

    pub fn edges(&self) -> [Segment; 3] {
        [
            Segment::new(self.a, self.b),
            Segment::new(self.b, self.c),
            Segment::new(self.c, self.a),
        ]
    }

    pub fn barycenter(&self) -> Point3<f64> {
        Point3::from((self.a.coords + self.b.coords + self.c.coords) / 3.0)
    }

    pub fn area(&self) -> f64 {
        let ab = self.b - self.a;
        let ac = self.c - self.a;
        0.5 * ab.cross(&ac).norm()
    }

    pub fn geogebra_string(&self) -> String {
        format!(
            "Polygon(({:}, {:}, {:}), ({:}, {:}, {:}), ({:}, {:}, {:}))",
            self.a.x,
            self.a.y,
            self.a.z,
            self.b.x,
            self.b.y,
            self.b.z,
            self.c.x,
            self.c.y,
            self.c.z
        )
    }
}

impl From<parry3d_f64::shape::Triangle> for Triangle {
    fn from(triangle: parry3d_f64::shape::Triangle) -> Self {
        Self {
            a: Point3::from(triangle.a),
            b: Point3::from(triangle.b),
            c: Point3::from(triangle.c),
        }
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use crate::primitive::segment::Segment;
    use crate::primitive::triangle::TriangleEdge;
    use crate::proptest::arbitrary_triangle;
    use proptest::proptest;

    proptest! {
        #[test]
        fn triangle_edge_returns_correct_segment(triangle in arbitrary_triangle()) {
            let edge_ab = triangle.edge(TriangleEdge::EdgeAB);
            let edge_bc = triangle.edge(TriangleEdge::EdgeBC);
            let edge_ca = triangle.edge(TriangleEdge::EdgeCA);
            assert_eq!(edge_ab, Segment::new(triangle.a, triangle.b));
            assert_eq!(edge_bc, Segment::new(triangle.b, triangle.c));
            assert_eq!(edge_ca, Segment::new(triangle.c, triangle.a));
        }
    }
}
