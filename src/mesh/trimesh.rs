use crate::primitive::aabb::Aabb;
use crate::primitive::triangle::{Triangle, TriangleEdge};
use nalgebra::Point3;
use ordered_float::OrderedFloat;
use std::collections::HashMap;
use union_find::{QuickUnionUf, UnionBySize, UnionFind};

/// Represents an index of a vertex in a mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VertexIndex(pub usize);

/// Represents an index of a triangle in a mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FaceIndex(pub usize);

/// Represents a triangle in a mesh, defined by three vertex indices.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct IndexTriangle {
    /// The first vertex index of the triangle.
    pub a: VertexIndex,
    /// The second vertex index of the triangle.
    pub b: VertexIndex,
    /// The third vertex index of the triangle.
    pub c: VertexIndex,
}

impl IndexTriangle {
    /// Creates a new `IndexTriangle` from three vertex indices.
    ///
    /// # Parameters
    /// - `a`: The first vertex index of the triangle.
    /// - `b`: The second vertex index of the triangle.
    /// - `c`: The third vertex index of the triangle.
    ///
    /// # Returns
    /// An `IndexTriangle` object representing the triangle.
    pub fn new(a: VertexIndex, b: VertexIndex, c: VertexIndex) -> Self {
        Self { a, b, c }
    }

    pub fn edges(&self) -> [[VertexIndex; 2]; 3] {
        [[self.a, self.b], [self.b, self.c], [self.c, self.a]]
    }

    pub fn edge_vertices(&self, triangle_edge: TriangleEdge) -> [VertexIndex; 2] {
        match triangle_edge {
            TriangleEdge::EdgeAB => [self.a, self.b],
            TriangleEdge::EdgeBC => [self.b, self.c],
            TriangleEdge::EdgeCA => [self.c, self.a],
        }
    }

    pub fn vertices(&self) -> [VertexIndex; 3] {
        [self.a, self.b, self.c]
    }
}

/// A point on the surface of a TriMesh.
#[derive(Clone, Copy, Debug)]
pub struct MeshSurfacePoint {
    /// The id of the face that the point lies on.
    pub face_id: FaceIndex,
    pub point: Point3<f64>,
}

/// Represents a triangle mesh, consisting of vertices and triangle indices.
#[derive(Debug, Clone)]
pub struct TriangleMesh {
    /// A list of vertices in the mesh.
    pub vertices: Vec<Point3<f64>>,
    /// A list of triangles, each defined by three vertex indices.
    pub indices: Vec<IndexTriangle>,
}

impl TriangleMesh {
    /// Creates a new `TriangleMesh` from a list of vertices and a list of triangles.
    ///
    /// # Parameters
    /// - `vertices`: A list of vertices in the mesh.
    /// - `indices`: A list of triangles, each defined by three vertex indices.
    ///
    /// # Returns
    /// A `TriangleMesh` object representing the mesh.
    pub fn new(vertices: Vec<Point3<f64>>, indices: Vec<IndexTriangle>) -> Self {
        Self { vertices, indices }
    }

    pub fn vertex(&self, vertex_index: VertexIndex) -> Point3<f64> {
        self.vertices[vertex_index.0]
    }

    pub fn triangle(&self, face_id: FaceIndex) -> IndexTriangle {
        self.indices[face_id.0]
    }

    pub fn euclidean_triangle(&self, face_id: FaceIndex) -> Triangle {
        let tri = self.triangle(face_id);
        Triangle {
            a: self.vertices[tri.a.0],
            b: self.vertices[tri.b.0],
            c: self.vertices[tri.c.0],
        }
    }

    pub fn triangles(&self) -> impl Iterator<Item = Triangle> + '_ {
        self.indices.iter().map(move |tri| Triangle {
            a: self.vertices[tri.a.0],
            b: self.vertices[tri.b.0],
            c: self.vertices[tri.c.0],
        })
    }

    pub fn surface_point_to_euclidean(&self, pt: MeshSurfacePoint) -> Point3<f64> {
        pt.point
    }

    #[cfg(test)]
    pub fn as_geogebra_string(&self) -> String {
        let mut result = String::new();
        result.push_str("{");

        for (i, tri) in self.triangles().enumerate() {
            if i > 0 {
                result.push_str(",");
            }
            result.push_str(&format!(
                "Polygon(({},{},{}),({},{},{}),({},{},{}))",
                tri.a.x, tri.a.y, tri.a.z, tri.b.x, tri.b.y, tri.b.z, tri.c.x, tri.c.y, tri.c.z
            ));
        }
        result.push_str("}");

        result
    }

    pub fn union(mut self, other: &Self) -> Self {
        let offset = self.vertices.len();

        self.vertices.extend(other.vertices.iter().cloned());
        self.indices
            .extend(other.indices.iter().map(|t| IndexTriangle {
                a: VertexIndex(t.a.0 + offset),
                b: VertexIndex(t.b.0 + offset),
                c: VertexIndex(t.c.0 + offset),
            }));

        self
    }

    pub fn rescaled(mut self, scale_factor: f64) -> Self {
        for vertex in &mut self.vertices {
            *vertex = Point3::new(
                vertex.x * scale_factor,
                vertex.y * scale_factor,
                vertex.z * scale_factor,
            );
        }
        self
    }
}

impl From<parry3d_f64::shape::TriMesh> for TriangleMesh {
    fn from(mesh: parry3d_f64::shape::TriMesh) -> Self {
        let vertices = mesh
            .vertices()
            .iter()
            .map(|v| Point3::new(v.x, v.y, v.z))
            .collect();
        let indices = mesh
            .indices()
            .iter()
            .map(|f| IndexTriangle {
                a: VertexIndex(f[0] as usize),
                b: VertexIndex(f[1] as usize),
                c: VertexIndex(f[2] as usize),
            })
            .collect();
        TriangleMesh::new(vertices, indices)
    }
}

impl From<&TriangleMesh> for parry3d_f64::shape::TriMesh {
    fn from(mesh: &TriangleMesh) -> Self {
        let vertices = mesh
            .vertices
            .iter()
            .map(|v| parry3d_f64::math::Point::new(v.x, v.y, v.z))
            .collect();
        let indices = mesh
            .indices
            .iter()
            .map(|f| [f.a.0 as u32, f.b.0 as u32, f.c.0 as u32])
            .collect();
        parry3d_f64::shape::TriMesh::new(vertices, indices).expect("Failed to create TriMesh")
    }
}

impl From<TriangleMesh> for parry3d_f64::shape::TriMesh {
    fn from(mesh: TriangleMesh) -> Self {
        Self::from(&mesh)
    }
}

/// Merge vertices that share the same position.
pub fn deduplicate_vertices(mesh: &TriangleMesh) -> TriangleMesh {
    let (old_to_new_index, new_vertices) = dedup_points_vec(&mesh.vertices);

    // Compute the new indices via the old-to-new index mapping.
    let new_indices: Vec<IndexTriangle> = mesh
        .indices
        .iter()
        .map(|tri| {
            IndexTriangle::new(
                VertexIndex(old_to_new_index[tri.a.0]),
                VertexIndex(old_to_new_index[tri.b.0]),
                VertexIndex(old_to_new_index[tri.c.0]),
            )
        })
        .collect();

    TriangleMesh::new(new_vertices, new_indices)
}

pub fn dedup_points_vec(points: &[Point3<f64>]) -> (Vec<usize>, Vec<Point3<f64>>) {
    // A hashmap that maps 3D points to their index in the new vertex list.
    let mut position_to_index = HashMap::new();
    // A list that maps old vertex indices to new vertex indices.
    let mut old_to_new_index = vec![0; points.len()];
    let mut new_vertices = Vec::new();

    for (i, vertex) in points.iter().enumerate() {
        // Annoyingly, we need to wrap the vertex in OrderedFloat to make it hashable.
        let key = [
            OrderedFloat(vertex.x),
            OrderedFloat(vertex.y),
            OrderedFloat(vertex.z),
        ];

        // If the vertex is already in the hashmap, use the existing index.
        let index = position_to_index.entry(key).or_insert_with(|| {
            // Otherwise, add the vertex to the new list and use the new index.
            let new_index = new_vertices.len();
            new_vertices.push(*vertex);
            new_index
        });

        // Update the mapping from old to new index.
        old_to_new_index[i] = *index;
    }

    // Make sure all the new indices are valid:
    for &index in &old_to_new_index {
        assert!(index < new_vertices.len());
    }

    (old_to_new_index, new_vertices)
}

pub struct ConnectedComponents {
    pub components: Vec<Vec<VertexIndex>>,
    pub triangle_components: Vec<Vec<FaceIndex>>,
}

impl ConnectedComponents {
    pub fn connected_components_of_mesh(mesh: &TriangleMesh) -> Self {
        // Create a union-find data structure to find connected components.
        let mut uf = QuickUnionUf::<UnionBySize>::new(mesh.vertices.len());

        // For each triangle, tell the union-find data structure that the vertices are connected.
        for i in 0..mesh.indices.len() {
            let tri = &mesh.indices[i];
            let a = tri.a.0;
            let b = tri.b.0;
            let c = tri.c.0;
            uf.union(a, b);
            uf.union(b, c);
            // We don't need to union a and c because we know that b is connected to both.
        }

        // Next, we will create a mapping from vertex index to connected component index.
        // We will re-index the connected components so that they are numbered from 0 to m-1,
        // where m is the number of connected components.

        let mut component_map = HashMap::new();
        let mut next_index = 0;

        let vertex_to_component = (0..mesh.vertices.len())
            .map(|i| {
                let component = uf.find(i);
                *component_map.entry(component).or_insert_with(|| {
                    let index = next_index;
                    next_index += 1;
                    index
                })
            })
            .collect::<Vec<_>>();

        let mut components = vec![Vec::new(); next_index];

        for (i, &component) in vertex_to_component.iter().enumerate() {
            components[component].push(VertexIndex(i));
        }

        let triangle_to_component: Vec<usize> = mesh
            .indices
            .iter()
            .map(|tri| {
                let a = vertex_to_component[tri.a.0];
                let b = vertex_to_component[tri.b.0];
                let c = vertex_to_component[tri.c.0];
                assert_eq!(a, b);
                assert_eq!(b, c);
                a
            })
            .collect();

        let mut triangle_components = vec![Vec::new(); next_index];

        for (i, &component) in triangle_to_component.iter().enumerate() {
            triangle_components[component].push(FaceIndex(i));
        }

        Self {
            components,
            triangle_components,
        }
    }
}

pub fn component_aabb(
    connected_components: &ConnectedComponents,
    mesh: &TriangleMesh,
    component_index: usize,
) -> Aabb {
    let vertices = connected_components.components[component_index]
        .iter()
        .map(|&v| mesh.vertex(v))
        .collect::<Vec<_>>();

    Aabb::around_points(vertices.iter().copied())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Extracts a specific connected component from a triangle mesh.
    ///
    /// # Parameters
    /// - `fruit_mesh`: The original triangle mesh from which the component is extracted.
    /// - `connected_components`: The connected components of the mesh, as computed by
    ///   `ConnectedComponents::connected_components_of_mesh`.
    /// - `ith_component`: The index of the connected component to extract.
    ///
    /// # Returns
    /// A new `TriangleMesh` containing only the vertices and triangles of the specified connected component.
    pub fn extract_connected_component(
        fruit_mesh: &TriangleMesh,
        connected_components: &ConnectedComponents,
        ith_component: usize,
    ) -> TriangleMesh {
        // An index translation map from the large mesh to the small mesh.
        let mut index_translation: HashMap<VertexIndex, VertexIndex> = HashMap::new();

        let mut vertices = Vec::new();

        // Map the vertices of the selected component to a new index and collect them.
        for vertex in &connected_components.components[ith_component] {
            let new_index = VertexIndex(vertices.len());
            vertices.push(fruit_mesh.vertex(*vertex));
            index_translation.insert(*vertex, new_index);
        }

        let mut indices = Vec::new();

        // Map the triangles of the selected component to the new vertex indices and collect them.
        for face in &connected_components.triangle_components[ith_component] {
            let tri = fruit_mesh.triangle(*face);
            let new_face = IndexTriangle::new(
                *index_translation.get(&tri.a).unwrap(),
                *index_translation.get(&tri.b).unwrap(),
                *index_translation.get(&tri.c).unwrap(),
            );
            indices.push(new_face);
        }

        // Return the new triangle mesh containing only the selected component.
        TriangleMesh { vertices, indices }
    }

    #[test]
    fn test_connected_components_extraction() {
        // Create three simple meshes at different locations

        // Mesh 1: Tetrahedron at origin
        let tetrahedron = TriangleMesh::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.5, 0.866, 0.0),
                Point3::new(0.5, 0.433, 0.816),
            ],
            vec![
                IndexTriangle::new(VertexIndex(0), VertexIndex(1), VertexIndex(2)),
                IndexTriangle::new(VertexIndex(0), VertexIndex(1), VertexIndex(3)),
                IndexTriangle::new(VertexIndex(1), VertexIndex(2), VertexIndex(3)),
                IndexTriangle::new(VertexIndex(0), VertexIndex(2), VertexIndex(3)),
            ],
        );

        // Mesh 2: Square (two triangles) at (10,0,0)
        let square = TriangleMesh::new(
            vec![
                Point3::new(10.0, 0.0, 0.0),
                Point3::new(11.0, 0.0, 0.0),
                Point3::new(11.0, 1.0, 0.0),
                Point3::new(10.0, 1.0, 0.0),
            ],
            vec![
                IndexTriangle::new(VertexIndex(0), VertexIndex(1), VertexIndex(2)),
                IndexTriangle::new(VertexIndex(0), VertexIndex(2), VertexIndex(3)),
            ],
        );

        // Mesh 3: Triangle at (0,0,10)
        let triangle = TriangleMesh::new(
            vec![
                Point3::new(0.0, 0.0, 10.0),
                Point3::new(1.0, 0.0, 10.0),
                Point3::new(0.5, 1.0, 10.0),
            ],
            vec![IndexTriangle::new(
                VertexIndex(0),
                VertexIndex(1),
                VertexIndex(2),
            )],
        );

        // Store the original meshes
        let original_meshes = vec![tetrahedron.clone(), square.clone(), triangle.clone()];
        let original_vertex_counts = original_meshes
            .iter()
            .map(|m| m.vertices.len())
            .collect::<Vec<_>>();
        let original_face_counts = original_meshes
            .iter()
            .map(|m| m.indices.len())
            .collect::<Vec<_>>();

        // Combine them into one mesh
        let combined_mesh = tetrahedron.union(&square).union(&triangle);

        // Verify the combined mesh has the expected number of vertices and triangles
        assert_eq!(combined_mesh.vertices.len(), 4 + 4 + 3);
        assert_eq!(combined_mesh.indices.len(), 4 + 2 + 1);

        // Find connected components
        let connected_components =
            ConnectedComponents::connected_components_of_mesh(&combined_mesh);

        // There should be 3 components
        assert_eq!(connected_components.components.len(), 3);

        // Extract the components
        let extracted_meshes: Vec<TriangleMesh> = (0..3)
            .map(|i| extract_connected_component(&combined_mesh, &connected_components, i))
            .collect();

        // Function to calculate the centroid of a mesh
        fn mesh_centroid(mesh: &TriangleMesh) -> Point3<f64> {
            let sum = mesh
                .vertices
                .iter()
                .fold(Point3::new(0.0, 0.0, 0.0), |acc, pt| {
                    Point3::new(acc.x + pt.x, acc.y + pt.y, acc.z + pt.z)
                });
            let count = mesh.vertices.len() as f64;
            Point3::new(sum.x / count, sum.y / count, sum.z / count)
        }

        // Match each extracted mesh with its original by comparing centroids
        for extracted in &extracted_meshes {
            let centroid = mesh_centroid(extracted);

            // Find which original mesh this corresponds to
            let mut found_match = false;
            for (i, original) in original_meshes.iter().enumerate() {
                let original_centroid = mesh_centroid(original);

                // Calculate distance between centroids
                let distance = (centroid - original_centroid).norm();

                if distance < 0.1 {
                    // Found a match
                    found_match = true;

                    // Verify counts
                    assert_eq!(extracted.vertices.len(), original_vertex_counts[i]);
                    assert_eq!(extracted.indices.len(), original_face_counts[i]);
                    break;
                }
            }

            assert!(
                found_match,
                "Could not match extracted mesh to any original"
            );
        }
    }
}
