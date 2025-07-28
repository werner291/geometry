use nalgebra::{Point3, Unit, Vector3};

#[derive(Debug, Copy, Clone)]
pub struct SurfacePoint {
    pub point: Point3<f64>,
    pub normal: Unit<Vector3<f64>>,
}
