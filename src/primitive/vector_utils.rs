use nalgebra::{Unit, UnitQuaternion, Vector3};

/// Compute an arbitrary 3D vector perpendicular to the given `Vector3<f64>`.
///
/// This function ensures numerical stability by comparing the absolute values of the components.
/// Credit: https://stackoverflow.com/a/11132720
///
/// # Parameters
/// - `vec`: A `Vector3<f64>` representing the input vector.
///
/// # Returns
/// A `Vector3<f64>` that is perpendicular to the input vector.
pub fn any_perpendicular(vec: Vector3<f64>) -> Vector3<f64> {
    if vec.z.abs() < vec.x.abs() {
        Vector3::new(vec.y, -vec.x, 0.0)
    } else {
        Vector3::new(0.0, -vec.z, vec.y)
    }
}

/// Rotates a vector to produce a second vector such that the angle between the two vectors
/// is `offset_angle` radians.
///
/// No guarantees are made about the direction of the resulting vector, except for
/// the provided offset angle.
///
/// However, this function takes an additional `plane_angle` parameter, which allows
/// sweeping the full range of possible vectors with the given offset angle.
///
/// # Parameters
///
/// - `original`: A `Unit<Vector3<f64>>` representing the original vector.
/// - `offset_angle`: The resulting angle in radians between the original vector and the resulting vector.
/// - `plane_angle`: An additional angle in radians to rotate the resulting vector around the original vector,
///                  providing a means to vary the resulting vector's direction.
pub fn offset_vector_by_angle(
    original: Unit<Vector3<f64>>,
    offset_angle: f64,
    plane_angle: f64,
) -> Unit<Vector3<f64>> {
    let plane_rotation = UnitQuaternion::from_axis_angle(&original, plane_angle);

    let offset_rotation = UnitQuaternion::from_axis_angle(
        // We don't really care about the direction of the perpendicular vector,
        // as long as it is perpendicular to the original vector.
        &Unit::new_normalize(any_perpendicular(*original)),
        offset_angle,
    );

    plane_rotation * offset_rotation * original
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::proptest::arbitrary_vector;
    use nalgebra::{Unit, Vector3};
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn any_perpendicular_is_perpendicular_to_input(vec in arbitrary_vector()) {
            if vec.norm() > 1e-10 {
                let result = any_perpendicular(vec);

                println!("Before: {}, after: {}", &vec, &result);

                println!("Dot is {}", result.dot(dbg!(&vec)).abs());
                prop_assert!(result.dot(&vec).abs() < 1e-10);
            }
        }

        #[test]
        fn offset_vector_by_angle_produces_correct_offset_angle(
            vec in arbitrary_vector(),
            offset_angle in 0.0..std::f64::consts::PI,
            plane_angle in 0.0..(2.0 * std::f64::consts::PI)
        ) {
            if vec.norm() > 1e-10 {
                let original = Unit::new_normalize(vec);
                let result = offset_vector_by_angle(original, offset_angle, plane_angle);
                let angle_between = result.dot(&original).acos();
                prop_assert!((angle_between - offset_angle).abs() < 1e-10);
            }
        }
    }
}
