use crate::mesh::trimesh::{FaceIndex, IndexTriangle, TriangleMesh, VertexIndex};
use crate::TriangleEdge;
use nalgebra::Point3;
use std::collections::HashMap;

/// Index type for half-edges.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct HalfEdgeIndex(pub usize);

/// Data associated with a vertex in the half-edge structure.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct VertexData {
    /// An incident half-edge originating from this vertex.
    pub incident_edge: HalfEdgeIndex,
}

/// Data associated with a half-edge in the half-edge structure.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HalfEdgeData {
    /// The vertex at the origin of this half-edge.
    pub origin_vertex: VertexIndex,
    /// The twin half-edge, which is the half-edge in the opposite direction.
    pub twin_edge: HalfEdgeIndex,
    /// The next half-edge in the face.
    pub next_edge: HalfEdgeIndex,
    /// The previous half-edge in the face.
    pub prev_edge: HalfEdgeIndex,
    /// The face to which this half-edge belongs.
    pub incident_face: FaceIndex,
}

/// Data associated with a face in the half-edge structure.
pub struct FaceData {
    /// An incident half-edge belonging to this face.
    pub incident_edge: HalfEdgeIndex,
}

/// A structure representing the half-edge topology of a mesh.
///
/// The half-edge data structure is a way to represent the connectivity of a mesh.
///
/// # Fields
///
/// * `vertices` - A vector of `VertexData` representing the vertices in the mesh.
/// * `half_edges` - A vector of `HalfEdgeData` representing the half-edges in the mesh.
/// * `faces` - A vector of `FaceData` representing the faces in the mesh.
pub struct HalfEdgeTopology {
    #[allow(dead_code)] // We'll be using this when treating the saddle points case.
    pub vertices: Vec<VertexData>,
    pub half_edges: Vec<HalfEdgeData>,
    pub faces: Vec<FaceData>,
}

impl HalfEdgeTopology {
    /// Constructs a `HalfEdgeTopology` from a given `TriangleMesh`.
    ///
    /// This function creates the half-edge data structure for the given triangle mesh by
    /// assigning indices to the half-edges, vertices, and faces, and establishing the
    /// connectivity between them.
    ///
    /// Precondition: the input mesh must be manifold.
    ///
    /// # Arguments
    ///
    /// * `mesh` - A reference to the `TriangleMesh` from which to construct the half-edge topology.
    ///
    /// # Returns
    ///
    /// A `HalfEdgeTopology` representing the connectivity of the given triangle mesh.
    ///
    /// # Panics
    ///
    /// This function will panic if the input mesh is not manifold.
    pub fn index_of_mesh(mesh: &TriangleMesh) -> Self {
        let index_assignment = assign_indices_to_halfedges(mesh);

        // Create a HashMap to store the incident edge for each vertex;
        // we use a hashmap for deduplication.
        let vertex_incident_edges: HashMap<VertexIndex, HalfEdgeIndex> = index_assignment
            .iter()
            .map(|([u, _], e)| (*u, *e))
            .collect();

        // For each vertex, look up the incident edge in the hashmap.
        let vertices = (0..mesh.vertices.len())
            .map(|i| VertexData {
                incident_edge: vertex_incident_edges[&VertexIndex(i)],
            })
            .collect();

        // For each face, create a FaceData object with the incident edge.
        // By convention, it will always be the AB edge of the triangle.
        let faces = mesh
            .indices
            .iter()
            .map(|tr| FaceData {
                incident_edge: index_assignment[&[tr.a, tr.b]],
            })
            .collect();

        // Collect the connectivity data for each half-edge.
        let half_edges = mesh
            .indices
            .iter()
            .enumerate()
            .flat_map(|(face_index, triangle)| {
                let edges = triangle.edges();
                let index_assignment = &index_assignment;
                // for each edge of the triangle...
                (0..3).map(move |i| HalfEdgeData {
                    // The origin vertex is the first vertex of the edge.
                    origin_vertex: edges[i][0],
                    // Get the twin by swapping the vertices and looking up the index.
                    twin_edge: index_assignment[&[edges[i][1], edges[i][0]]],
                    // Get the next one by looking up the next edge in the triangle.
                    next_edge: index_assignment[&[edges[(i + 1) % 3][0], edges[(i + 1) % 3][1]]],
                    // Previous edge is the next of the next edge since a triangle has 3 edges.
                    prev_edge: index_assignment[&[edges[(i + 2) % 3][0], edges[(i + 2) % 3][1]]],
                    incident_face: FaceIndex(face_index),
                })
            })
            .collect();

        Self {
            vertices,
            half_edges,
            faces,
        }
    }

    /// Returns the next half-edge in the face.
    pub fn next(&self, e: HalfEdgeIndex) -> HalfEdgeIndex {
        self.half_edges[e.0].next_edge
    }

    /// Returns the twin half-edge.
    pub fn twin(&self, e: HalfEdgeIndex) -> HalfEdgeIndex {
        self.half_edges[e.0].twin_edge
    }

    /// Returns the vertex at the origin of the half-edge.
    pub fn origin(&self, e: HalfEdgeIndex) -> VertexIndex {
        self.half_edges[e.0].origin_vertex
    }

    /// Returns the face to which the half-edge belongs.
    pub fn incident_face(&self, e: HalfEdgeIndex) -> FaceIndex {
        self.half_edges[e.0].incident_face
    }

    /// Returns the incident half-edge of a face.
    pub fn incident_edge_of_face(&self, f: FaceIndex) -> HalfEdgeIndex {
        self.faces[f.0].incident_edge
    }

    /// Returns an iterator over the half-edges of a face.
    pub fn face_halfedges(&self, f: FaceIndex) -> impl Iterator<Item = HalfEdgeIndex> + '_ {
        let first_edge = self.incident_edge_of_face(f);
        std::iter::successors(Some(first_edge), move |&e| {
            if self.next(e) == first_edge {
                None
            } else {
                Some(self.next(e))
            }
        })
    }

    /// Returns the vertices of an edge.
    pub fn edge_vertices(&self, e: HalfEdgeIndex) -> [VertexIndex; 2] {
        [self.origin(e), self.origin(self.twin(e))]
    }

    /// Returns the points of an edge.
    pub fn edge_points(&self, e: HalfEdgeIndex, mesh: &TriangleMesh) -> [Point3<f64>; 2] {
        let [u, v] = self.edge_vertices(e);
        [mesh.vertex(u), mesh.vertex(v)]
    }

    /// Returns an iterator over the half-edges incident to a given vertex.
    ///
    /// This function starts from an incident half-edge of the vertex and follows
    /// the half-edges around the vertex in a counter-clockwise direction until it
    /// returns to the starting half-edge.
    ///
    /// # Arguments
    ///
    /// * `v` - The index of the vertex for which to find the incident half-edges.
    ///
    /// # Returns
    ///
    /// An iterator over the `HalfEdgeIndex` values of the half-edges incident to the vertex.
    pub fn halfedges_incident_to_vertex(
        &self,
        v: VertexIndex,
    ) -> impl Iterator<Item = HalfEdgeIndex> + '_ {
        // The first edge is just whichever incident edge the vertex has.
        let first_edge = self.vertices[v.0].incident_edge;

        std::iter::successors(Some(first_edge), move |&e| {
            // By taking the twin then the next, we can follow the edges around the vertex.
            let next = self.next(self.twin(e));
            if next == first_edge {
                None
            } else {
                Some(next)
            }
        })
    }
}

/// Assigns indices to the half-edges of a given `TriangleMesh`.
///
/// This function iterates over the edges of the triangles in the mesh and assigns a unique index
/// to each half-edge. The indices are stored in a `HashMap` where the key is a pair of vertex indices
/// representing the half-edge, and the value is the assigned `HalfEdgeIndex`.
///
/// # Arguments
///
/// * `mesh` - A reference to the `TriangleMesh` for which to assign half-edge indices.
///
/// # Returns
///
/// A `HashMap` where the key is a pair of `VertexIndex` representing the half-edge, and the value is the assigned `HalfEdgeIndex`.
fn assign_indices_to_halfedges(mesh: &TriangleMesh) -> HashMap<[VertexIndex; 2], HalfEdgeIndex> {
    // Create a HashMap to store the index assignment for each half-edge.
    let index_assignment: HashMap<[VertexIndex; 2], HalfEdgeIndex> = mesh
        // Iterate over the triangles in the mesh.
        .indices
        .iter()
        // Flatten the edges of all triangles in the mesh into an iterator of edges.
        .flat_map(IndexTriangle::edges)
        // Add an index.
        .enumerate()
        // Flip it, so the key is the pair of vertices and the value is the index.
        .map(|(i, v)| (v, HalfEdgeIndex(i)))
        // Collect the results into a HashMap.
        .collect();

    // Ensure that the number of assigned indices matches the expected number of half-edges.
    assert_eq!(index_assignment.len(), mesh.indices.len() * 3);

    // Return the index assignment.
    index_assignment
}

/// Convert a half-edge to a `TriangleEdge` enum value representing one of the edges of a triangle.
/// That is, we're trying to identify the edge of the triangle that the half-edge corresponds to.
///
/// # Arguments
///
/// * `half_edge` - The index of the half-edge in the topology.
/// * `tri_mesh` - The triangle mesh containing the vertices.
/// * `topology` - The half-edge topology of the triangle mesh.
///
/// # Returns
///
/// A `TriangleEdge` enum value representing the edge of the triangle that the half-edge corresponds to.
///
/// # Panics
///
/// This function will panic if the half-edge does not match any of the triangle edges.
pub fn halfedge_to_ab_bc_ca(
    half_edge: HalfEdgeIndex,
    tri_mesh: &TriangleMesh,
    topology: &HalfEdgeTopology,
) -> TriangleEdge {
    let [v0, v1] = topology.edge_vertices(half_edge);
    let indices = tri_mesh.triangle(topology.incident_face(half_edge));

    if indices.a == v0 && indices.b == v1 {
        TriangleEdge::EdgeAB
    } else if indices.b == v0 && indices.c == v1 {
        TriangleEdge::EdgeBC
    } else if indices.c == v0 && indices.a == v1 {
        TriangleEdge::EdgeCA
    } else {
        panic!("Half-edge does not match any of the triangle edges.")
    }
}

/// Converts a `TriangleEdge` to a `HalfEdgeIndex`.
///
/// This function identifies the half-edge in the half-edge topology that corresponds to the given
/// edge of the triangle.
///
/// # Arguments
///
/// * `edge` - The `TriangleEdge` enum value representing one of the edges of the triangle.
/// * `face` - The index of the face in the topology.
/// * `tri_mesh` - The triangle mesh containing the vertices.
/// * `topology` - The half-edge topology of the triangle mesh.
///
/// # Returns
///
/// A `HalfEdgeIndex` representing the half-edge in the topology that corresponds to the given edge.
///
/// # Panics
///
/// This function will panic if the edge origin vertex does not match any of the face vertices.
pub fn ab_bc_ca_to_halfedge(
    edge: TriangleEdge,
    face: FaceIndex,
    tri_mesh: &TriangleMesh,
    topology: &HalfEdgeTopology,
) -> HalfEdgeIndex {
    let edge_origin_vertex = match edge {
        TriangleEdge::EdgeAB => tri_mesh.triangle(face).a,
        TriangleEdge::EdgeBC => tri_mesh.triangle(face).b,
        TriangleEdge::EdgeCA => tri_mesh.triangle(face).c,
    };

    let he0 = topology.incident_edge_of_face(face);
    let he1 = topology.next(he0);
    let he2 = topology.next(he1);

    if topology.origin(he0) == edge_origin_vertex {
        he0
    } else if topology.origin(he1) == edge_origin_vertex {
        he1
    } else if topology.origin(he2) == edge_origin_vertex {
        he2
    } else {
        panic!("Edge origin vertex does not match any of the face vertices.")
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::proptest::arbitrary_convex_hull;
    use proptest::proptest;

    proptest! {

        #[test]
        fn test_halfedge_to_ab_bc_ca_cycle_is_identity(chull in arbitrary_convex_hull(1000)) {
            let topology = HalfEdgeTopology::index_of_mesh(&chull);

            for face in 0..chull.indices.len() {
                let face = FaceIndex(face);
                for half_edge in topology.face_halfedges(face) {
                    let edge = halfedge_to_ab_bc_ca(half_edge, &chull, &topology);
                    let half_edge_back = ab_bc_ca_to_halfedge(edge, face, &chull, &topology);
                    assert_eq!(
                        half_edge, half_edge_back,
                        "Cycle conversion failed for half_edge: {:?}",
                        half_edge
                    );
                }
            }
        }

        #[test]
        fn convex_hull_topology_is_valid(chull in arbitrary_convex_hull(1000)) {
            let topology = HalfEdgeTopology::index_of_mesh(&chull);

            for half_edge in 0..topology.half_edges.len() {
                let half_edge = HalfEdgeIndex(half_edge);
                let twin_of_twin = topology.twin(topology.twin(half_edge));
                assert_eq!(half_edge, twin_of_twin);
                let triple_next = topology.next(topology.next(topology.next(half_edge)));
                assert_eq!(half_edge, triple_next);
            }

            for face in 0..topology.faces.len() {
                let face = FaceIndex(face);
                let incident_edge = topology.incident_edge_of_face(face);
                let face_edges: Vec<_> = topology.face_halfedges(face).collect();
                assert!(
                    face_edges.contains(&incident_edge),
                    "Incident edge of face validity failed for face: {:?}",
                    face
                );

                // The AB half edge should be the one that the face data points to.
                let ab_vertices = topology.edge_vertices(incident_edge);
                let face_vertices = chull.triangle(face).edges();
                assert_eq!(ab_vertices, [face_vertices[0][0], face_vertices[0][1]]);
            }

            for vertex in 0..topology.vertices.len() {
                // The source vertex of the incident edge should be the vertex itself.
                let vertex = VertexIndex(vertex);
                let incident_edge = topology.vertices[vertex.0].incident_edge;
                assert_eq!(topology.origin(incident_edge), vertex);

                for half_edge in topology.halfedges_incident_to_vertex(vertex) {
                    assert_eq!(topology.origin(half_edge), vertex);
                }
            }
        }
    }
}
