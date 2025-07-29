//! This module contains a Rust implementation of the shortest path algorithm on triangular meshes. (without saddle points!)
//!
//! It is based on the documentation https://doc.cgal.org/5.6.2/Surface_mesh_shortest_path/index.html#title19

#![allow(unused)]

use crate::intersection::{parametric_line_intersection, RaySegmentIntersection};
use crate::intersection::{ray_segment_intersection, EdgePoint};
use crate::mesh::halfedge::{
    ab_bc_ca_to_halfedge, halfedge_to_ab_bc_ca, HalfEdgeIndex, HalfEdgeTopology,
};
use crate::mesh::trimesh::MeshSurfacePoint;
use crate::mesh::trimesh::{TriangleMesh, VertexIndex};
use crate::mesh_shortest_paths::PropagationResult::{SingleWindow, SplitWindows};
use crate::primitive::infinite_triangle::InfiniteTriangle;
use crate::primitive::ray::Ray;
use crate::primitive::segment::Segment;
// Note: These imports are commented out as they reference modules from the main project
// that are not part of the standalone geometry crate
// use crate::shell_path_planning::from_point_and_normal::StateFromPointAndNormal;
// use crate::visualizations::common::render_utils::ConvertThreeD;
// use crate::visualizations::common::render_utils::{Visualize, line_segments_mesh};
use itertools::Itertools;
use nalgebra::{distance, Const, OPoint, OVector, Point3, Unit, Vector3};
use ordered_float::OrderedFloat;
use parry3d_f64::shape::TriMesh;
use petgraph::prelude::UnGraph;
use std::cmp::PartialEq;
use std::collections::BinaryHeap;

/// A datastructure that forms a tree-shaped roadmap of the surface of a TriMesh;
/// it is similar to a Dijkstra shortest-path tree, and can be used to quickly
/// find geodesics along the surface of a mesh, starting from some source point.
pub struct SequenceTree {
    /// A vector of windows; the root windows are always at the start of this vector
    pub windows: Vec<WindowFaceTraversal>,
    /// A mapping of, for each face, which visibility windows cross into it.
    pub face_window_mapping: Vec<Vec<usize>>,
}

pub struct GeodesicSegment<'a> {
    pub traversal: &'a WindowFaceTraversal,
    pub entry_point: Option<Point3<f64>>,
    pub exit_point: Option<Point3<f64>>,
}

impl SequenceTree {
    /// Retrace a sequence of visibility windows from a given point.
    ///
    /// This method looks up which the visibility windows in the face that the point lies on,
    /// finds the window through which the predicted geodesic distance would be shortest,
    /// and then retraces the visibility windows back to the root.
    ///
    /// This method returns an iterator over the visibility windows, starting from the leaf
    /// window and ending at the root window.
    ///
    /// # Arguments
    ///
    /// * `tm` - The triangle mesh that the point lies on.
    /// * `q` - The point to retrace from.
    ///
    /// # Returns
    ///
    /// An iterator over the visibility windows, starting from the leaf window and ending at the root window,
    /// or `None` if the point is not in a visibility window.
    pub fn retrace(
        &self,
        tm: &TriangleMesh,
        q: &MeshSurfacePoint,
    ) -> Option<impl Iterator<Item = GeodesicSegment>> {
        {
            // Sanity check: q must lie in the plane of the face.
            let tri = tm.euclidean_triangle(q.face_id);
            let normal = tri.unit_normal();
            let delta = q.point - tri.a;
            let dot = delta.dot(&normal);
            assert!(dot.abs() < 1.0e-6);
        }

        // Look up the visibility windows this face is associated with:
        // Find the best window:
        let windows = &self.face_window_mapping[q.face_id.0];

        let window_index = (windows)
            .iter()
            .filter(|w| {
                let coverage = &self.windows[**w].coverage;
                let pt = &tm.surface_point_to_euclidean(*q);
                let delta = pt - coverage.origin;

                // Check if the point is in the plane of the triangle.
                assert!((coverage.normal().dot(&delta)).abs() < 1.0e-5);

                let n1 = coverage.left_dir.cross(&coverage.normal());
                let n2 = coverage.normal().cross(&coverage.right_dir);

                const MARGIN: f64 = 1.0e-6;

                delta.dot(&n1) <= MARGIN && delta.dot(&n2) <= MARGIN
            })
            .min_by_key(|w| OrderedFloat(self.windows[**w].distance_in_window(q, tm)))
            .expect("The point should be in a visibility window.");
        // TODO ?;

        let mut lat_point_euclidean = q.point;

        // Retrace the path:
        Some(
            std::iter::successors(Some(*window_index), |w| self.windows[*w].parent).map(move |w| {
                let is_first = w == *window_index;
                let traversal = &self.windows[w];

                // lat_point_euclidean must lie in the plane of the coverage window:
                let d = traversal
                    .coverage
                    .normal()
                    .dot(&(lat_point_euclidean - traversal.coverage.origin));
                assert!(d.abs() < 1.0e-6);

                let ray = Ray::new(
                    traversal.coverage.origin,
                    Unit::new_normalize(lat_point_euclidean - traversal.coverage.origin),
                );

                let exit_point = (!is_first).then(|| {
                    let exit_segment = Segment {
                        a: traversal.face_exit.left_point,
                        b: traversal.face_exit.right_point,
                    };

                    let [t1, t2] = ray_segment_intersection(ray, exit_segment)
                        .intersection()
                        .expect("The ray should intersect the exit segment.");
                    ray.at(t1)
                });

                let entry_point = traversal.parent.map(|parent| {
                    let parent = &self.windows[parent];

                    let entry_segment = Segment {
                        a: parent.face_exit.left_point,
                        b: parent.face_exit.right_point,
                    };

                    // a and b must be in the plane of the parent window:
                    {
                        let d1 = parent
                            .coverage
                            .normal()
                            .dot(&(entry_segment.a - parent.coverage.origin));
                        let d2 = parent
                            .coverage
                            .normal()
                            .dot(&(entry_segment.b - parent.coverage.origin));
                        assert!(d1.abs() < 1.0e-6);
                        assert!(d2.abs() < 1.0e-6);
                    }

                    let [t1, t2] = ray_segment_intersection(ray, entry_segment)
                        .intersection()
                        .expect("The ray should intersect the entry segment.");

                    let entry = ray.at(t1);
                    let pt2 = ray.at(t2);

                    // Entry must be on the edge of the parent window:
                    {
                        let edge = Segment {
                            a: parent.face_exit.left_point,
                            b: parent.face_exit.right_point,
                        };
                        let d = edge.distance_to_point(entry);
                        let d2 = edge.distance_to_point(pt2);
                        assert!(d < 1.0e-6);
                    }

                    // must lie in the plane of the parent window:
                    {
                        let d = parent
                            .coverage
                            .normal()
                            .dot(&(entry - parent.coverage.origin));
                        assert!(d.abs() < 1.0e-6);
                    }

                    entry
                });

                if let Some(entry_point) = entry_point {
                    lat_point_euclidean = entry_point;
                }

                GeodesicSegment {
                    traversal,
                    entry_point,
                    exit_point,
                }
            }),
        )
    }

    pub fn extract_quad(&self, traversal_id: usize) -> [Point3<f64>; 4] {
        let traversal = &self.windows[traversal_id];
        if let Some(parent) = traversal.parent {
            let parent = &self.windows[parent];
            [
                parent.face_exit.left_point,
                parent.face_exit.right_point,
                traversal.face_exit.left_point,
                traversal.face_exit.right_point,
            ]
        } else {
            [
                traversal.coverage.origin,
                traversal.coverage.origin,
                traversal.face_exit.right_point,
                traversal.face_exit.left_point,
            ]
        }
    }
}

/// A struct representing the portion of an edge that a visibility window exits through.
///
/// This struct contains information about the half-edge through which the face exits,
/// the intersection parameters for the left and right bounding rays, and the points
/// of intersection in 3D space.
#[derive(Debug)]
pub struct FaceExit {
    /// The index of the half-edge in the topology.
    pub half_edge: HalfEdgeIndex,
    /// The intersection parameter for the left bounding ray.
    left_t: f64,
    /// The intersection parameter for the right bounding ray.
    right_t: f64,
    /// The point of intersection for the left bounding ray in 3D space.
    pub left_point: Point3<f64>,
    /// The point of intersection for the right bounding ray in 3D space.
    pub right_point: Point3<f64>,
}

/// A struct representing a visibility window traversing a face, exiting through an edge.
///
/// Visibility windows are ordered by the length of the shortest geodesic starting within them.
#[derive(Debug)]
pub struct WindowFaceTraversal {
    /// The infinite triangle that this window covers.
    pub coverage: InfiniteTriangle,
    /// The face entry from some parent window, if any.
    pub parent: Option<usize>,
    /// The face exit to some child window.
    pub face_exit: FaceExit,
    /// The depth of the tree branch.
    tree_depth: usize,
}

impl WindowFaceTraversal {
    fn new(
        coverage: InfiniteTriangle,
        parent: Option<usize>,
        face_exit: FaceExit,
        tree_depth: usize,
    ) -> Self {
        Self {
            coverage,
            parent,
            face_exit,
            tree_depth,
        }
    }

    /// The length of the geodesic starting from the given point.
    fn distance_in_window(&self, pt: &MeshSurfacePoint, mesh: &TriangleMesh) -> f64 {
        // Check that this point is for the right face.
        distance(&self.coverage.origin, &mesh.surface_point_to_euclidean(*pt))
    }

    /// The length of the shortest geodesic between the source of the window,
    /// and any point on the edge that the window exits through.
    fn heuristic_distance(&self) -> f64 {
        let left = self.face_exit.left_point;
        let right = self.face_exit.right_point;
        let origin = self.coverage.origin;

        Segment::new(left, right).distance_to_point(origin)
    }
}

/// Computes a vector perpendicular to the given half-edge in the plane of the triangle,
/// pointing towards the vertex of the triangle opposite the edge.
///
/// # Arguments
///
/// * `half_edge` - The index of the half-edge in the topology.
/// * `tri_mesh` - The triangle mesh containing the vertices.
/// * `topology` - The half-edge topology of the triangle mesh.
///
/// # Returns
///
/// A `Vector3<f64>` representing the perpendicular vector.
fn halfedge_perpendicular(
    half_edge: HalfEdgeIndex,
    tri_mesh: &TriangleMesh,
    topology: &HalfEdgeTopology,
) -> Vector3<f64> {
    // Find the other two half-edges:
    let e2 = topology.next(half_edge);
    let e3 = topology.next(e2);

    // Look up the vertices of the half-edges.
    let v1 = topology.origin(half_edge);
    let v2 = topology.origin(e2);
    let v3 = topology.origin(e3);

    // And the corresponding points. Given how half-edge topology is constructed,
    // v1 and v2 are the vertices of the half-edge, and v3 is the vertex opposite the half-edge.
    let v1 = tri_mesh.vertex(v1);
    let v2 = tri_mesh.vertex(v2);
    let v3 = tri_mesh.vertex(v3);

    // Compute a vector parallel to the edge, normalized.
    let basis_x = (v2 - v1).normalize();

    // Project v3 onto the edge.
    let v3_proj = v1 + (v3 - v1).dot(&basis_x) * basis_x;

    // Compute a vector perpendicular to the edge, pointing towards v3.

    (v3 - v3_proj).normalize()
}

/// A struct representing an expansion head in the priority queue, referring
/// to a leaf of the sequence tree that is a candidate for further expansion/propagation.
///
/// This struct contains the ID of the window and the heuristic distance
/// used for ordering in the priority queue.
#[derive(Debug)]
pub struct ExpansionHead {
    pub window_id: usize,
    pub heuristic_distance: f64,
}

impl PartialEq for ExpansionHead {
    /// Checks if two `ExpansionHead` instances are equal based on their heuristic distance.
    fn eq(&self, other: &Self) -> bool {
        self.heuristic_distance == other.heuristic_distance
    }
}

impl Eq for ExpansionHead {}

impl PartialOrd for ExpansionHead {
    /// Compares two `ExpansionHead` instances based on their heuristic distance.
    /// The ordering is reversed because the priority queue is a max-heap, so a lower
    /// heuristic distance is considered "greater" or "higher priority".
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ExpansionHead {
    /// Compares two `ExpansionHead` instances based on their heuristic distance.
    /// The ordering is reversed because the priority queue is a max-heap, so a lower
    /// heuristic distance is considered "greater" or "higher priority".
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .heuristic_distance
            .partial_cmp(&self.heuristic_distance)
            .expect("Heuristic distance should always be finite.")
    }
}

pub struct MeshShortestPathsHooks<'a> {
    pub on_root_window_created: &'a dyn Fn(usize, &WindowFaceTraversal),
    pub on_visit_head: &'a dyn Fn(&ExpansionHead, &WindowFaceTraversal),
    pub head_discarded_max_depth: &'a dyn Fn(&ExpansionHead, &WindowFaceTraversal),
    pub head_discarded_xin_wang: &'a dyn Fn(&ExpansionHead, &WindowFaceTraversal),
    pub xin_wang_evaluated: &'a dyn Fn(&WindowFaceTraversal, XinWangResult),
    pub window_propagated: &'a dyn Fn(&ExpansionHead, &WindowFaceTraversal, &PropagationResult),
    pub vertex_best_distance_updated: &'a dyn Fn(VertexIndex, f64, f64, &WindowFaceTraversal),
}

/// Build a sequence tree for a given source point on a mesh.
///
/// A sequence tree is an acceleration structure that allows for the quick computation of geodesics
/// on a triangular mesh that start or end at the given source point.
///
/// Note that this makes this a one-to-many shortest path algorithm, so if you are computing many
/// geodesics, it is beneficial to batch them together based on shared endpoints.
///
/// A sequence tree is analogous to the shortest-path tree in Dijkstra's algorithm, but for
/// the surface of a triangular mesh, rather than a graph or a grid.
///
/// # Arguments
///
/// * `source_point` -  The point on the mesh to build the sequence tree from.
/// * `tri_mesh` -      The triangle mesh to build the sequence tree on.
/// * `topology` -      The half-edge topology of the triangle mesh.
///
/// # Returns
///
/// A `SequenceTree` representing the sequence tree.
pub fn build_sequence_tree(
    source_point: MeshSurfacePoint,
    tri_mesh: &TriangleMesh,
    topology: &HalfEdgeTopology,
    hooks: Option<MeshShortestPathsHooks>,
) -> SequenceTree {
    /////////////////////////
    // Datastructure init: //
    /////////////////////////

    // A vector of visibility windows.
    // Invariant:
    // - `windows` is ordered such that the parent of a window is always before the child, if the child has been instantiated.
    // - This implies that the root windows are always at the start of the vector.
    // - This vector is always appended to; indices are never invalidated but window children may be added.
    let mut traversals = Vec::new();

    // For every vertex in the mesh, the best distance to that vertex from the source point.
    let mut vertex_best_distance = init_vertex_distances(tri_mesh);

    // A lookup vector of vectors of visibility windows, one per face.
    let mut windows_in_face = init_window_faces(tri_mesh);

    // Priority queue of candidate windows.
    let mut expansion_heads: BinaryHeap<ExpansionHead> = BinaryHeap::new();

    //////////////////////////////
    // Geometry initialization: //
    //////////////////////////////

    let at_vertex = tri_mesh.indices[source_point.face_id.0]
        .vertices()
        .into_iter()
        .find(|&v| distance(&tri_mesh.vertex(v), &source_point.point) < 1.0e-6);

    if let Some(at_vertex) = at_vertex {
        // This is a special case; we want to create a window for opposite edge of the vertex.
        let vertex = tri_mesh.vertex(at_vertex);

        for he in topology.halfedges_incident_to_vertex(at_vertex) {
            let he = topology.next(he);
            create_root_window_through_edge(
                source_point.point,
                tri_mesh,
                topology,
                &hooks,
                &mut traversals,
                &mut windows_in_face,
                &mut expansion_heads,
                he,
            );
        }
    } else {
        // We might be on one of the edges of the face:
        let edge = topology.face_halfedges(source_point.face_id).find(|&he| {
            let [v1, v2] = topology.edge_points(he, tri_mesh);
            let edge = Segment::new(v1, v2);

            let distance = edge.distance_to_point(source_point.point);

            distance < 1.0e-5
        });

        if let Some(edge) = edge {
            // The windows will exit through the two edges opposite this edge,
            // and through those opposite the twin of this edge.

            let twin = topology.twin(edge);

            let through_edges = [
                topology.next(edge),
                topology.next(topology.next(edge)),
                topology.next(twin),
                topology.next(topology.next(twin)),
            ];

            for he in through_edges {
                create_root_window_through_edge(
                    source_point.point,
                    tri_mesh,
                    topology,
                    &hooks,
                    &mut traversals,
                    &mut windows_in_face,
                    &mut expansion_heads,
                    he,
                );
            }
        } else {
            // Create a window for each edge of the starting face:
            for he in topology.face_halfedges(source_point.face_id) {
                create_root_window_through_edge(
                    source_point.point,
                    tri_mesh,
                    topology,
                    &hooks,
                    &mut traversals,
                    &mut windows_in_face,
                    &mut expansion_heads,
                    he,
                );
            }
        }
    }

    ////////////////////////
    // Propagation phase: //
    ////////////////////////

    // The maximum depth of the sequence tree is the number of faces in the mesh.
    let max_depth = tri_mesh.indices.len();

    while let Some(head) = expansion_heads.pop() {
        let window = &traversals[head.window_id];

        if let Some(hooks) = hooks.as_ref() {
            (hooks.on_visit_head)(&head, window)
        }

        // First, check if we must discard the window.
        if window.parent.is_some() {
            if window.tree_depth >= max_depth {
                if let Some(hooks) = hooks.as_ref() {
                    (hooks.head_discarded_max_depth)(&head, window)
                }
                continue;
            }

            if xin_wang_rule(
                tri_mesh,
                &vertex_best_distance,
                topology.twin(
                    traversals[window
                        .parent
                        .expect("Earlier is_some check should ensure a parent here")]
                    .face_exit
                    .half_edge,
                ),
                window,
                topology,
                &hooks,
            ) {
                if let Some(hooks) = hooks.as_ref() {
                    (hooks.head_discarded_xin_wang)(&head, window)
                }
                continue;
            }
        }

        // Enqueue candidate children windows for propagation.
        let children = propagate(head.window_id, window, tri_mesh, topology);

        if let Some(hooks) = hooks.as_ref() {
            (hooks.window_propagated)(&head, window, &children)
        }

        match children {
            PropagationResult::SingleWindow(w) => {
                validate_left_right(tri_mesh, topology, &w);
                enqueue_window(
                    tri_mesh,
                    topology,
                    &mut traversals,
                    &mut windows_in_face,
                    &mut expansion_heads,
                    w,
                );
                // TODO: try applying the Xin-Wang rule here as well to prune early.
            }
            PropagationResult::SplitWindows {
                left,
                splitting_vertex,
                right,
            } => {
                validate_left_right(tri_mesh, topology, &left);
                validate_left_right(tri_mesh, topology, &right);

                // Update the best distance to the splitting vertex.
                let best_distance = vertex_best_distance[splitting_vertex.0];
                assert_eq!(left.coverage.origin, right.coverage.origin);
                let new_distance =
                    distance(&tri_mesh.vertex(splitting_vertex), &left.coverage.origin);

                if new_distance < best_distance {
                    if let Some(hooks) = hooks.as_ref() {
                        (hooks.vertex_best_distance_updated)(
                            splitting_vertex,
                            new_distance,
                            best_distance,
                            window,
                        )
                    }
                    vertex_best_distance[splitting_vertex.0] = new_distance;
                }

                // TODO: is this needed if the distance is longer?
                for w in [left, right] {
                    enqueue_window(
                        tri_mesh,
                        topology,
                        &mut traversals,
                        &mut windows_in_face,
                        &mut expansion_heads,
                        w,
                    );
                }
            }
        }
    }

    // Return the sequence tree.
    SequenceTree {
        windows: traversals,
        face_window_mapping: windows_in_face,
    }
}

fn enqueue_window(
    mesh: &TriangleMesh,
    topology: &HalfEdgeTopology,
    traversals: &mut Vec<WindowFaceTraversal>,
    windows_in_face: &mut Vec<Vec<usize>>,
    expansion_heads: &mut BinaryHeap<ExpansionHead>,
    w: WindowFaceTraversal,
) {
    validate_left_right(mesh, topology, &w);

    let window_id = traversals.len();
    expansion_heads.push(ExpansionHead {
        window_id,
        heuristic_distance: w.heuristic_distance(),
    });
    // Add the window to the face it is on.

    {
        // Normals must match:
        let face_normal = mesh
            .euclidean_triangle(topology.incident_face(w.face_exit.half_edge))
            .unit_normal();
        let window_normal = w.coverage.normal().normalize();
        let dot = face_normal.dot(&window_normal);
        assert!(dot > 0.9999);
    }

    windows_in_face[topology.incident_face(w.face_exit.half_edge).0].push(window_id);
    traversals.push(w);
}

fn validate_left_right(mesh: &TriangleMesh, topology: &HalfEdgeTopology, w: &WindowFaceTraversal) {
    // Quick validation step:
    let exit_edge_points = topology.edge_points(w.face_exit.half_edge, mesh);

    // Their direction vectors should be the same.
    let exit_dir = (exit_edge_points[1] - exit_edge_points[0]).normalize();
    let exit_segment_dir = (w.face_exit.right_point - w.face_exit.left_point).normalize();
    let dot = exit_dir.dot(&exit_segment_dir);
    assert!(dot > 0.9999);
}

fn create_root_window_through_edge(
    source_point: Point3<f64>,
    tri_mesh: &TriangleMesh,
    topology: &HalfEdgeTopology,
    hooks: &Option<MeshShortestPathsHooks>,
    traversals: &mut Vec<WindowFaceTraversal>,
    windows_in_face: &mut Vec<Vec<usize>>,
    expansion_heads: &mut BinaryHeap<ExpansionHead>,
    he: HalfEdgeIndex,
) {
    // Create a visibility window through the current half-edge.
    let window_id = traversals.len();
    let [v1, v2] = topology.edge_points(he, tri_mesh);
    let left_delta = v1 - source_point;
    let right_delta = v2 - source_point;

    // Check that they're not 0-norm vectors.
    assert!(left_delta.norm_squared() > 1.0e-6);
    assert!(right_delta.norm_squared() > 1.0e-6);
    // TODO: There's an edge case here where the norm *is* zero.
    // Note to later self: it's important that we cover the entire span of directions.
    let window = WindowFaceTraversal {
        coverage: InfiniteTriangle::new(
            source_point,
            Unit::new_normalize(left_delta),
            Unit::new_normalize(right_delta),
        ),
        parent: None,
        face_exit: FaceExit {
            half_edge: he,
            left_t: 0.0,
            right_t: 1.0,
            left_point: v1,
            right_point: v2,
        },
        tree_depth: 0,
    };

    {
        let left_ray = window.coverage.left_ray();
        let right_ray = window.coverage.right_ray();
        let edge = Segment::new(v1, v2);

        ray_segment_intersection(left_ray, edge)
            .intersection()
            .expect("The left ray should intersect the edge.");
        ray_segment_intersection(right_ray, edge)
            .intersection()
            .expect("The right ray should intersect the edge.");

        // Same for the twin of that edge:
        let twin = topology.twin(he);
        let twin_edge = Segment::new(v2, v1);
        ray_segment_intersection(left_ray, twin_edge)
            .intersection()
            .expect("The left ray should intersect the twin edge.");
        ray_segment_intersection(right_ray, twin_edge)
            .intersection()
            .expect("The right ray should intersect the twin edge.");
    }

    // Associate the window with the face it is on.
    windows_in_face[topology.incident_face(he).0].push(window_id);

    {
        // Normals must match:
        let face_normal = tri_mesh
            .euclidean_triangle(topology.incident_face(he))
            .unit_normal();
        let window_normal = window.coverage.normal().normalize();
        let dot = face_normal.dot(&window_normal);
        assert!(dot > 0.9999);
    }

    // Enqueue the window as a candidate for expansion.
    expansion_heads.push(ExpansionHead {
        window_id,
        heuristic_distance: window.heuristic_distance(),
    });

    if let Some(hooks) = hooks.as_ref() {
        (hooks.on_root_window_created)(window_id, &window)
    }

    validate_left_right(tri_mesh, topology, &window);

    // Add the created window to the list of windows.
    traversals.push(window);
}

#[derive(Debug)]
pub enum PropagationResult {
    SingleWindow(WindowFaceTraversal),
    SplitWindows {
        left: WindowFaceTraversal,
        splitting_vertex: VertexIndex,
        right: WindowFaceTraversal,
    },
}

impl PropagationResult {
    pub fn windows(&self) -> Vec<&WindowFaceTraversal> {
        match self {
            PropagationResult::SingleWindow(w) => vec![w],
            PropagationResult::SplitWindows { left, right, .. } => vec![left, right],
        }
    }
}

/// Rotates a visibility window over an edge to the neighboring triangle.
///
/// This function computes the new position and directions of a visibility window
/// after it crosses an edge to the neighboring triangle.
///
/// # Arguments
///
/// * `parent_window` - The visibility window in the originating triangle.
/// * `he` - The index of the half-edge in the topology.
/// * `tri_mesh` - The triangle mesh containing the vertices.
/// * `topology` - The half-edge topology of the triangle mesh.
///
/// # Returns
///
/// An `InfiniteTriangle` representing the rotated visibility window in the neighboring triangle.
fn rotate_window_over_edge(
    parent_window: &InfiniteTriangle,
    he: HalfEdgeIndex,
    tri_mesh: &TriangleMesh,
    topology: &HalfEdgeTopology,
) -> InfiniteTriangle {
    let perp = halfedge_perpendicular(he, tri_mesh, topology);
    let twin = topology.twin(he);
    let twin_perp = halfedge_perpendicular(twin, tri_mesh, topology);
    let left_dir = bend_vector(perp, twin_perp, *parent_window.left_dir);
    let right_dir = bend_vector(perp, twin_perp, *parent_window.right_dir);
    let origin = bend_point(
        perp,
        twin_perp,
        topology.edge_points(he, tri_mesh),
        parent_window.origin,
    );

    // Make sure the origin is at the same distance from the edge as the parent window.
    assert!({
        let edge = topology.edge_points(he, tri_mesh);
        let delta_1 = origin - edge[0];
        let delta_2 = parent_window.origin - edge[0];
        ((delta_1.norm() - delta_2.norm()).abs() < 1.0e-6)
    });

    let new_triangle_normal = tri_mesh
        .euclidean_triangle(topology.incident_face(twin))
        .unit_normal();

    let ld = new_triangle_normal.dot(&left_dir);

    // Make sure the left and right directions lie in the plane of the triangle.
    assert!(ld.abs() < 1.0e-6);

    let rd = new_triangle_normal.dot(&right_dir);

    // Make sure the left and right directions lie in the plane of the triangle.
    assert!(rd.abs() < 1.0e-6);

    let tri = InfiniteTriangle {
        origin,
        left_dir: Unit::new_normalize(left_dir),
        right_dir: Unit::new_normalize(right_dir),
    };

    // The rays should intersect the edge:
    #[cfg(debug_assertions)]
    {
        // The face normal and new_coverage normals should align:
        assert!(new_triangle_normal.dot(&tri.normal().normalize()) > 0.9999);

        let edge_points = topology.edge_points(he, tri_mesh);
        let edge = Segment::new(edge_points[0], edge_points[1]);
        let left_ray = tri.left_ray();
        let right_ray = tri.right_ray();

        let ray_segment_intersection_tolerance: f64 = 1.0e-6;
        let left_intersection =
            parametric_line_intersection(left_ray.origin, *left_ray.dir, edge.a, edge.b - edge.a);
        let right_intersection =
            parametric_line_intersection(right_ray.origin, *right_ray.dir, edge.a, edge.b - edge.a);

        if let Some([t1, t2]) = left_intersection {
            assert!((-ray_segment_intersection_tolerance
                ..=1.0 + ray_segment_intersection_tolerance)
                .contains(&t2));
        }
        if let Some([t1, t2]) = right_intersection {
            assert!((-ray_segment_intersection_tolerance
                ..=1.0 + ray_segment_intersection_tolerance)
                .contains(&t2));
        }

        // If the above are None they're just parallel. That's fine; this is just a santiy check.
    }

    tri
}

const FLOAT_ASSERTION_MARGIN: f64 = 1.0e-10;

/// Propagate a given visibility window through whatever triangle it's traversing,
/// producing between one and two child windows depending on whether it occupies
/// the vertex opposite the edge it enters through.
///
/// Note: this function does not check whether the resulting windows are to be pruned,
/// either according to the Xin-Wang rule, tree depth, or one-vertex-one-split rule.
///
/// See:
///
/// ![img](src/geometry/ray_triangle_exit_cases.svg)
///
/// # Arguments
///
/// * `parent_id` - The ID of the parent window being propagated.
/// * `parent_window` - The parent window being propagated.
/// * `tri_mesh` - The triangle mesh that we're working with.
/// * `topology` - The half-edge topology of the triangle mesh.
///
/// # Returns
///
/// A `PropagationResult` representing the result of the propagation, containing either
/// one or two windows.
fn propagate(
    parent_id: usize,
    parent_window: &WindowFaceTraversal,
    tri_mesh: &TriangleMesh,
    topology: &HalfEdgeTopology,
) -> PropagationResult {
    // The half-edge through which the parent window exits its face.
    let he = parent_window.face_exit.half_edge;

    // First, double-check that the rays of the parent coverage window intersect the edge
    // that the parent window exits through. This *should* be the case, but it's good to be
    // sure as we're dealing with a semi-complex algorithm here.
    debug_assert!({
        let edge_points = topology.edge_points(he, tri_mesh);
        let edge = Segment::new(edge_points[0], edge_points[1]);
        !matches!(
            ray_segment_intersection(parent_window.coverage.left_ray(), edge),
            RaySegmentIntersection::NoIntersection
        ) && !matches!(
            ray_segment_intersection(parent_window.coverage.right_ray(), edge),
            RaySegmentIntersection::NoIntersection
        )
    });

    // Look up the twin of the parent's exit edge, which will be the entry edge for the next face.
    let twin = topology.twin(he);

    // Apply a rotation of the parent's coverage window over the edge to align it with the neighboring face.
    let mut new_coverage = rotate_window_over_edge(&parent_window.coverage, he, tri_mesh, topology);
    // Look up the index of the face we're entering.
    let new_face = topology.incident_face(twin);

    // The normal of the face we just entered, and the normal of the new coverage window should align.
    debug_assert!(
        tri_mesh
            .euclidean_triangle(new_face)
            .unit_normal()
            .dot(&new_coverage.normal().normalize())
            .abs()
            > 0.9999
    );

    // Extract the bounding rays of the new coverage window.
    let left_ray = new_coverage.left_ray();
    let right_ray = new_coverage.right_ray();

    // Assert that the rays are perpendicular to the triangle normal.
    debug_assert!({
        let new_face_normal = tri_mesh.euclidean_triangle(new_face).unit_normal();
        new_face_normal.dot(&left_ray.dir).abs() < 1.0e-6
            && new_face_normal.dot(&right_ray.dir).abs() < 1.0e-6
    });

    // Now, we're going to define the "left" and "right" edges of the new face.
    //
    // These are the two edges other than the edge through which we entered the face;
    // we may be exiting through either of them, or both of them.

    // The left edge is the next edge of the entry edge.
    let left_edge = topology.next(twin);
    let left_edge_points = topology.edge_points(left_edge, tri_mesh);
    let left_edge_segment = Segment::new(left_edge_points[0], left_edge_points[1]);

    // The right edge is the next of the left edge; they thus share a vertex.
    let right_edge = topology.next(left_edge);
    let right_edge_points = topology.edge_points(right_edge, tri_mesh);
    let right_edge_segment = Segment::new(right_edge_points[0], right_edge_points[1]);

    let shared_vertex = topology.origin(right_edge);
    let shared_point = tri_mesh.vertex(shared_vertex);

    // println!("LSeg = {}", left_edge_segment.geogebra_string());
    // println!("RSeg = {}", right_edge_segment.geogebra_string());
    // println!("LRay = {}", left_ray.geogebra_string());
    // println!("RRay = {}", right_ray.geogebra_string());
    // println!("Tri = {}", tri_mesh.euclidean_triangle(new_face).geogebra_string());

    let delta = &shared_point - new_coverage.origin;

    // Check if the point is in the plane of the triangle.
    debug_assert!((new_coverage.normal().dot(&delta)).abs() < 1.0e-5);

    let cov_normal = new_coverage.left_dir.cross(&new_coverage.right_dir);
    let n_left = new_coverage.left_dir.cross(&delta);
    let n_right = delta.cross(&new_coverage.right_dir);

    let left_ray_through_left_edge = (n_left.dot(&cov_normal) > 0.0).then(|| {
        parametric_line_intersection(
            left_ray.origin,
            *left_ray.dir,
            left_edge_segment.a,
            left_edge_segment.b - left_edge_segment.a,
        )
        .map(|t| t[1])
    });

    let right_ray_through_right_edge = (n_right.dot(&cov_normal) > 0.0).then(|| {
        parametric_line_intersection(
            right_ray.origin,
            *right_ray.dir,
            right_edge_segment.a,
            right_edge_segment.b - right_edge_segment.a,
        )
        .map(|t| t[1])
    });

    match (left_ray_through_left_edge, right_ray_through_right_edge) {
        (Some(left), Some(right)) => {
            // The left ray passes on the left of the shared vertex, the right one on the right.
            // Thus, the shared vertex splits the window in two
            SplitWindows {
                left: through_edge_with_ts(
                    tri_mesh,
                    topology,
                    parent_id,
                    parent_window,
                    &mut new_coverage,
                    left_edge,
                    [left.unwrap_or(0.0), 1.0],
                ),
                splitting_vertex: shared_vertex,
                right: through_edge_with_ts(
                    tri_mesh,
                    topology,
                    parent_id,
                    parent_window,
                    &mut new_coverage,
                    right_edge,
                    [0.0, right.unwrap_or(1.0)],
                ),
            }
        }
        (Some(left), None) => {
            // The left passes through the left, the right also passes on the left.
            let right_ray_left_intersection = parametric_line_intersection(
                right_ray.origin,
                *right_ray.dir,
                left_edge_segment.a,
                left_edge_segment.b - left_edge_segment.a,
            )
            .filter(|t| t[0] >= 0.0)
            .map(|t| t[1])
            .unwrap_or(1.0);

            SingleWindow(through_edge_with_ts(
                tri_mesh,
                topology,
                parent_id,
                parent_window,
                &mut new_coverage,
                left_edge,
                [
                    left.expect("The left ray should (cleanly) intersect the left edge."),
                    right_ray_left_intersection,
                ],
            ))
        }
        (None, Some(right)) => {
            // The left passes on the right, the right also passes on the right.
            let left_ray_right_intersection = parametric_line_intersection(
                left_ray.origin,
                *left_ray.dir,
                right_edge_segment.a,
                right_edge_segment.b - right_edge_segment.a,
            )
            .filter(|t| t[0] >= 0.0)
            .map(|t| t[1])
            .unwrap_or(1.0);

            SingleWindow(through_edge_with_ts(
                tri_mesh,
                topology,
                parent_id,
                parent_window,
                &mut new_coverage,
                right_edge,
                [
                    left_ray_right_intersection,
                    right.expect("The right ray should (cleanly) intersect the right edge."),
                ],
            ))
        }
        (None, None) => panic!(
            "This should never occur, as the shared vertex is simultanoesly left of the left ray and right of the right ray."
        ),
    }
}

fn through_edge_with_ts(
    tri_mesh: &TriangleMesh,
    topology: &HalfEdgeTopology,
    parent_id: usize,
    parent_window: &WindowFaceTraversal,
    new_coverage: &mut InfiniteTriangle,
    edge: HalfEdgeIndex,
    ts: [f64; 2],
) -> WindowFaceTraversal {
    // look up the edge
    let edge_points = topology.edge_points(edge, tri_mesh);

    assert!(ts[0] >= -0.01 && ts[0] <= 1.01);
    assert!(ts[1] >= -0.01 && ts[1] <= 1.01);

    let left_point = edge_points[0] + ts[0].clamp(0.0, 1.0) * (edge_points[1] - edge_points[0]);
    let right_point = edge_points[0] + ts[1].clamp(0.0, 1.0) * (edge_points[1] - edge_points[0]);

    let w = WindowFaceTraversal {
        coverage: InfiniteTriangle {
            origin: new_coverage.origin,
            left_dir: Unit::new_normalize(left_point - new_coverage.origin),
            right_dir: Unit::new_normalize(right_point - new_coverage.origin),
        },
        parent: Some(parent_id),
        face_exit: FaceExit {
            half_edge: edge,
            left_t: ts[0],
            right_t: ts[1],
            left_point,
            right_point,
        },
        tree_depth: parent_window.tree_depth + 1,
    };

    validate_left_right(tri_mesh, topology, &w);

    w
}

fn heuristic_distance_from_exit(new_coverage: &InfiniteTriangle, face_exit: &FaceExit) -> f64 {
    let ld = (face_exit.left_point - new_coverage.origin).norm();
    let rd = (face_exit.right_point - new_coverage.origin).norm();
    ld.min(rd)
}

fn init_window_faces(tri_mesh: &TriangleMesh) -> Vec<Vec<usize>> {
    tri_mesh
        .indices
        .iter()
        .map(|_| Vec::new())
        .collect::<Vec<Vec<usize>>>()
}

fn init_vertex_distances(tri_mesh: &TriangleMesh) -> Vec<f64> {
    tri_mesh
        // Once for each vertex in the mesh.
        .vertices
        .iter()
        // Initialize to infinity, to indicate that no distance has been found yet.
        .map(|_| f64::INFINITY)
        .collect::<Vec<f64>>()
}

/// Bends a point around an edge by rotating it about a given edge.
///
/// The tangents are assumed to be normalized, perpendicular to the edge,
/// pointing towards their respective opposite vertex.
///
/// # Arguments
///
/// * `tangent_1` - The tangent vector of the originating triangle (see `halfedge_perpendicular`).
/// * `tangent_2` - The tangent vector of the neighboring triangle (see `halfedge_perpendicular`).
/// * `edge_vertices` - The vertices of the edge shared by the two triangles.
/// * `pt` - The point to be bent.
///
/// # Returns
///
/// A new point representing the bent point.
fn bend_point(
    tangent_1: Vector3<f64>,
    tangent_2: Vector3<f64>,
    edge_vertices: [Point3<f64>; 2],
    pt: Point3<f64>,
) -> Point3<f64> {
    edge_vertices[0] + bend_vector(tangent_1, tangent_2, pt - edge_vertices[0])
}

/// Computes the new direction of a vector after bending it around an edge,
/// given a vector tangent to the originating triangle and a vector tangent
/// to the neighboring triangle, perpendicular to the edge.
///
/// The tangents are assumed to be normalized, and perpendicular to the edge,
/// each pointing towards their respective opposite vertex.
///
/// # Arguments
///
/// * `tangent_1` - The tangent vector of the originating triangle (see `halfedge_perpendicular`).
/// * `tangent_2` - The tangent vector of the neighboring triangle (see `halfedge_perpendicular`).
/// * `parent_left` - The vector to be bent.
///
/// # Returns
///
/// A new `Vector3<f64>` representing the bent vector.
fn bend_vector(
    tangent_1: Vector3<f64>,
    tangent_2: Vector3<f64>,
    parent_left: Vector3<f64>,
) -> Vector3<f64> {
    // Project the vector onto the first tangent.
    let l_factor = parent_left.dot(&tangent_1);
    // Subtract the tangent (scaled by the projection factor) from the vector, then add the tangent of the other triangle.
    // (we actually subtract since the two tangents point in opposite directions).
    parent_left - l_factor * tangent_1 - l_factor * tangent_2
}

#[derive(Debug)]
pub struct XinWangResult {
    pub comparisons_made: [(Point3<f64>, VertexIndex, bool); 3],
}

/// Applies the Xin-Wang rule (Theorem 3.2) to determine if a visibility window should be discarded.
///
/// THEOREM 3.2. Let w be a window that enter triangle v1 v2 v3 through
/// edge v1 v3. Assume that d1 , d2 , d3 are respectively the minimum-so-
/// far distance at the three vertices v 1 , v2 , v 3 . Then w can’t define a
/// shortest sequence if
/// d + |I B| > d1 + |v1 B|,
/// or
/// d + |I A| > d2 + |v2 A|,
/// or
/// d + |I A| > d3 + |v3 A|.
///
/// # Arguments
///
/// * `tri_mesh` - The triangle mesh containing the vertices.
/// * `vertex_best_distance` - A vector of the best known distances to each vertex in the mesh.
/// * `entry_edge` - The edge through which the window enters the triangle; only applies if not a root window.
/// * `window` - The visibility window being evaluated.
///
/// # Returns
///
/// `true` if the visibility window should be discarded, `false` otherwise.
fn xin_wang_rule(
    tri_mesh: &TriangleMesh,
    vertex_best_distance: &[f64],
    entry_edge: HalfEdgeIndex,
    window: &WindowFaceTraversal,
    topo: &HalfEdgeTopology,
    hooks: &Option<MeshShortestPathsHooks>,
) -> bool {
    // Here comes the tricky bit, because the Xin-Wang paper only describes one half of a symmetric case,
    // and we're going to have to fill in the other half.

    let exit_edge = window.face_exit.half_edge;

    let [right_entry, left_entry] = topo.edge_vertices(entry_edge);
    let [left_exit, right_exit] = topo.edge_vertices(exit_edge);

    // |I A| + d (d is always 0 because we don't consider saddle vertices)
    let dist_exit_left = (window.coverage.origin - window.face_exit.left_point).norm();
    // |I B| + d
    let dist_exit_right = (window.coverage.origin - window.face_exit.right_point).norm();

    if left_entry == left_exit {
        // The left vertex is shared; there are two vertices on the right side.
        // This means that the distances we want to compare are:
        //      B to v1, and A to v2 and v3.
        let i_left = left_entry;
        let i_right_1 = right_entry;
        let i_right_2 = right_exit;

        let best_left = vertex_best_distance[i_left.0];
        let best_right1 = vertex_best_distance[i_right_1.0];
        let best_right2 = vertex_best_distance[i_right_2.0];

        let left: Point3<f64> = tri_mesh.vertex(i_left);
        let right_1 = tri_mesh.vertex(i_right_1);
        let right_2 = tri_mesh.vertex(i_right_2);

        // Quick check, do we have left and right mixed up?
        assert!(
            (left - window.face_exit.left_point).norm()
                <= (left - window.face_exit.right_point).norm()
        );
        assert!(
            (right_2 - window.face_exit.right_point).norm()
                <= (right_2 - window.face_exit.left_point).norm()
        );

        let dist_left = (left - window.face_exit.right_point).norm() + best_left;
        let dist_right1 = (right_1 - window.face_exit.left_point).norm() + best_right1;
        let dist_right2 = (right_2 - window.face_exit.left_point).norm() + best_right2;

        let discard = dist_exit_right > dist_left
            || dist_exit_left > dist_right1
            || dist_exit_left > dist_right2;

        if let Some(hooks) = hooks.as_ref() {
            (hooks.xin_wang_evaluated)(
                window,
                XinWangResult {
                    comparisons_made: [
                        (
                            window.face_exit.right_point,
                            i_left,
                            dist_exit_right > dist_left,
                        ),
                        (
                            window.face_exit.left_point,
                            i_right_1,
                            dist_exit_left > dist_right1,
                        ),
                        (
                            window.face_exit.left_point,
                            i_right_2,
                            dist_exit_left > dist_right2,
                        ),
                    ],
                },
            )
        }

        discard
    } else {
        // The right vertex is shared; there are two vertices on the left side.
        // This means that the distances we want to compare are:
        //      A to v1, and B to v2 and v3.
        let left_1 = left_entry;
        let left_2 = left_exit;
        let right = right_entry;

        let best_left1 = vertex_best_distance[left_1.0];
        let best_left2 = vertex_best_distance[left_2.0];
        let best_right = vertex_best_distance[right.0];

        let vertex_left_1 = tri_mesh.vertex(left_1);
        let vertex_left_2 = tri_mesh.vertex(left_2);
        let vertex_right = tri_mesh.vertex(right);

        // Quick check, do we have left and right mixed up?
        assert!(
            (vertex_left_2 - window.face_exit.left_point).norm()
                <= (vertex_left_2 - window.face_exit.right_point).norm()
        );
        assert!(
            (vertex_right - window.face_exit.right_point).norm()
                <= (vertex_right - window.face_exit.left_point).norm()
        );

        let dist_left1 = (vertex_left_1 - window.face_exit.right_point).norm() + best_left1;
        let dist_left2 = (vertex_left_2 - window.face_exit.right_point).norm() + best_left2;
        let dist_right = (vertex_right - window.face_exit.left_point).norm() + best_right;

        let discard = dist_exit_right > dist_left1
            || dist_exit_right > dist_left2
            || dist_exit_left > dist_right;

        if let Some(hooks) = hooks.as_ref() {
            (hooks.xin_wang_evaluated)(
                window,
                XinWangResult {
                    comparisons_made: [
                        (
                            window.face_exit.right_point,
                            left_1,
                            dist_exit_right > dist_left1,
                        ),
                        (
                            window.face_exit.right_point,
                            left_2,
                            dist_exit_right > dist_left2,
                        ),
                        (
                            window.face_exit.left_point,
                            right,
                            dist_exit_left > dist_right,
                        ),
                    ],
                },
            )
        }

        discard
    }
}

fn retrace_shell_path(
    chull: &TriangleMesh,
    topology: &HalfEdgeTopology,
    start: MeshSurfacePoint,
    end: &MeshSurfacePoint,
    st: &SequenceTree,
) -> Vec<MeshSurfacePoint> {
    let mut path = st
        .retrace(chull, &end)
        .expect("The convex hull should be connected.")
        .flat_map(|traversal| {
            let face_id = topology.incident_face(traversal.traversal.face_exit.half_edge);

            let entry = traversal
                .entry_point
                .unwrap_or(chull.surface_point_to_euclidean(start));
            let exit = traversal
                .exit_point
                .unwrap_or(chull.surface_point_to_euclidean(*end));

            [
                MeshSurfacePoint {
                    face_id,
                    point: exit, // Since we're retracing, gotta put the exit point first.
                },
                MeshSurfacePoint {
                    face_id,
                    point: entry,
                },
            ]
        })
        .collect::<Vec<_>>();

    path.reverse();

    path
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::mesh::convex_hull::compute_chull;
    use crate::mesh::trimesh::FaceIndex;
    use crate::primitive::triangle::TriangleEdge;
    use crate::proptest::{
        arbitrary_convex_hull, arbitrary_point_on_arbitrary_chull,
        arbitrary_point_pair_on_arbitrary_chull,
    };
    use petgraph::prelude::UnGraph;
    use proptest::{prop_assert, proptest};
    use rand::prelude::StdRng;
    use rand::{Rng, SeedableRng};

    pub fn triangle_from_halfedge(
        p0: HalfEdgeIndex,
        p1: &TriangleMesh,
        p2: &HalfEdgeTopology,
    ) -> crate::primitive::triangle::Triangle {
        p1.euclidean_triangle(p2.incident_face(p0))
    }

    enum EdgeGraphPoint {
        Vertex(VertexIndex),
        Start,
        Goal,
    }

    type EdgeGraph = UnGraph<EdgeGraphPoint, f64>;

    fn shortest_path_via_edges(
        starting_point: MeshSurfacePoint,
        ending_point: MeshSurfacePoint,
        chull: &TriangleMesh,
    ) -> f64 {
        let mut graph = EdgeGraph::new_undirected();

        let mesh_vertices = chull
            .vertices
            .iter()
            .enumerate()
            .map(|(i, _)| graph.add_node(EdgeGraphPoint::Vertex(VertexIndex(i))))
            .collect::<Vec<_>>();

        for triangle in chull.indices.iter() {
            for edge in triangle.edges() {
                let p1 = chull.vertex(edge[0]);
                let p2 = chull.vertex(edge[1]);
                let dist = distance(&p1, &p2);
                graph.add_edge(mesh_vertices[edge[0].0], mesh_vertices[edge[1].0], dist);
            }
        }

        // Add the start and end points:
        let start = graph.add_node(EdgeGraphPoint::Start);
        let end = graph.add_node(EdgeGraphPoint::Goal);

        // Connect the start and end points to the mesh by creating an edge to every vertex in the face:
        for vertex in chull.triangle(starting_point.face_id).vertices() {
            graph.add_edge(
                start,
                mesh_vertices[vertex.0],
                distance(&starting_point.point, &chull.vertex(vertex)),
            );
        }

        for vertex in chull.triangle(ending_point.face_id).vertices() {
            graph.add_edge(
                end,
                mesh_vertices[vertex.0],
                distance(&ending_point.point, &chull.vertex(vertex)),
            );
        }

        // Run Dijkstra's algorithm to find the shortest path between the start and end points.
        let result = petgraph::algo::dijkstra(&graph, start, Some(end), |e| *e.weight());

        result[&end]
    }

    proptest! {

        #[test]
        fn test_halfedge_perpendicular_tangent(chull in arbitrary_convex_hull(1000)) {
            let topology = HalfEdgeTopology::index_of_mesh(&chull);

            for he in 0..topology.half_edges.len() {
                let he = HalfEdgeIndex(he);
                let triangle = triangle_from_halfedge(he, &chull, &topology);

                let perp = halfedge_perpendicular(he, &chull, &topology);

                // It must be normalized.
                assert!((perp - perp.normalize()).norm() < 1.0e-8);

                // Must be perpendicular to the triangle normal, i.e. in the tangent plane.
                let normal = triangle.unit_normal();
                assert!(perp.dot(&normal).abs() < 1.0e-8);

                // Must be perpendicular to the edge:
                let edge = topology.edge_points(he, &chull);
                let edge_dir = (edge[1] - edge[0]).normalize();
                assert!(perp.dot(&edge_dir).abs() < 1.0e-8);

                // Must point towards the opposite vertex:
                let v_opposite = chull.vertex(topology.origin(topology.next(topology.next(he))));
                assert!(perp.dot(&(v_opposite - edge[0])) > 0.0);
            }
        }

        #[ignore] // <- Test is outdated since I added more asserstions in the method.
        #[test]
        fn test_rotate_window_over_edge_twin(chull in arbitrary_convex_hull(1000)) {
            let topology = HalfEdgeTopology::index_of_mesh(&chull);
            let mut rng = StdRng::seed_from_u64(42);

            for he in 0..topology.half_edges.len() {
                let he = HalfEdgeIndex(he);
                let twin = topology.twin(he);

                let triangle = triangle_from_halfedge(he, &chull, &topology);
                let normal = triangle.unit_normal();

                // We don't really care about the window being a valid traversal, we only want it to be tangent to the triangle.
                let origin = Point3::new(
                    rng.random_range(-10.0..=10.0),
                    rng.random_range(-10.0..=10.0),
                    rng.random_range(-10.0..=10.0),
                );
                let left_dir = Vector3::new(
                    rng.random_range(-1.0..=1.0),
                    rng.random_range(-1.0..=1.0),
                    rng.random_range(-1.0..=1.0),
                );
                let right_dir = Vector3::new(
                    rng.random_range(-1.0..=1.0),
                    rng.random_range(-1.0..=1.0),
                    rng.random_range(-1.0..=1.0),
                );
                // Project all on the triangle plane:
                let origin = origin - (origin - triangle.a).dot(&normal) * *normal;
                assert!((origin - triangle.a).dot(&normal).abs() < 1e-8);
                let left_dir = Unit::new_normalize(left_dir - left_dir.dot(&normal) * *normal);
                let right_dir = Unit::new_normalize(right_dir - right_dir.dot(&normal) * *normal);

                let window = InfiniteTriangle {
                    origin,
                    left_dir,
                    right_dir,
                };

                let rotated_once = rotate_window_over_edge(&window, he, &chull, &topology);

                // Check that the resulting window lies in the tangent plane of the twin triangle.
                let new_triangle = triangle_from_halfedge(twin, &chull, &topology);
                let new_normal = new_triangle.unit_normal();
                let d_from_t = (rotated_once.origin - new_triangle.a)
                    .dot(&new_normal)
                    .abs();
                assert!(rotated_once.left_dir.dot(&new_normal).abs() < 1e-8);
                assert!(rotated_once.right_dir.dot(&new_normal).abs() < 1e-8);
                assert!(d_from_t < 1e-8);

                // Check that the projection of the origin on the edge is the same as the original origin.
                // Comparing the dot products should be enough.
                let edge = topology.edge_points(he, &chull);
                let d1 = (edge[1] - edge[0]).normalize().dot(&(rotated_once.origin - edge[0]));
                let d2 = (edge[1] - edge[0]).normalize().dot(&(window.origin - edge[0]));
                prop_assert!( (d1 - d2).abs() < 1e-8 );
                // Same for the directions:
                prop_assert!(
                    ((edge[1] - edge[0]).dot(&rotated_once.left_dir)
                        - (edge[1] - edge[0]).dot(&window.left_dir))
                    .abs()
                        < 1e-8
                );
                prop_assert!(
                    ((edge[1] - edge[0]).dot(&rotated_once.right_dir)
                        - (edge[1] - edge[0]).dot(&window.right_dir))
                    .abs()
                        < 1e-8
                );

                let rotated_twice = rotate_window_over_edge(&rotated_once, twin, &chull, &topology);

                // Check that the window is the same as the original after rotating twice
                prop_assert!((rotated_twice.origin - window.origin).norm() < 1e-8);
                prop_assert!((*rotated_twice.left_dir - *window.left_dir).norm() < 1e-8);
                prop_assert!((*rotated_twice.right_dir - *window.right_dir).norm() < 1e-8);
            }
        }

        #[test]
        fn all_windows_touched_by_tree((source, chull) in arbitrary_point_on_arbitrary_chull(1000)) {
            println!("Chull has {} vertices and {} faces.", chull.vertices.len(), chull.indices.len());
            println!("Chull = {}", chull.as_geogebra_string());
            println!("Source = ({}, {}, {})", source.point.x, source.point.y, source.point.z);

            let topology = HalfEdgeTopology::index_of_mesh(&chull);

            let st = build_sequence_tree(
                source,
                &chull,
                &topology,
                None
            );

            // For all facets:
            for f in 0..chull.indices.len() {
                // it must have at least one window in the face:
                let face_windows = &st.face_window_mapping[f];

                println!("Window triangle: {:?}", chull.triangle(FaceIndex(f)));
                println!("Points: {}", chull.euclidean_triangle(FaceIndex(f)).geogebra_string());
                println!("Windows in face: {:?}", face_windows);

                if face_windows.is_empty() {
                    println!("No windows in face.");
                }

                prop_assert!(!face_windows.is_empty());

                // All windows' normal must be aligned with the face normal.
                let face_normal = chull.euclidean_triangle(FaceIndex(f)).unit_normal();
                for window in face_windows.iter() {

                    // It must backpoint to the face.
                    assert_eq!(topology.incident_face(st.windows[*window].face_exit.half_edge), FaceIndex(f));

                    let normal = st.windows[*window].coverage.normal().normalize();
                    let dt = normal.dot(&face_normal);
                    assert!(dt.abs() > 0.9999);
                }
            }
        }

        #[test]
        fn test_convex_hull_versus_edge_graph(([source_point, goal_point], chull) in arbitrary_point_pair_on_arbitrary_chull(1000)) {

            println!("Testing source point: {:?}", source_point);
            println!("Testing goal point: {:?}", goal_point);
            println!("Chull has {} vertices and {} faces.", chull.vertices.len(), chull.indices.len());
            println!("Chull = {}", chull.as_geogebra_string());

            let edge_graph_distance = shortest_path_via_edges(source_point, goal_point, &chull);

            println!("Edge graph distance: {}", edge_graph_distance);

            // This is the "big test".
            // Specifically, we compare the length of paths found by the sequence tree to the length
            // of paths found by just running Dijkstra's algorithm on the edge graph.

            let topology = HalfEdgeTopology::index_of_mesh(&chull);

            let st = build_sequence_tree(
                source_point,
                &chull,
                &topology,
                None
            );

            println!("Start=({}, {}, {})", source_point.point.x, source_point.point.y, source_point.point.z);
            println!("Goal=({}, {}, {})", goal_point.point.x, goal_point.point.y, goal_point.point.z);

            let path = retrace_shell_path(&chull, &topology, source_point, &goal_point, &st);
            let path_length = path.windows(2).map(|w| distance(&w[0].point, &w[1].point)).sum::<f64>();

            println!("Path has {} points.", path.len());
            println!("Path length: {}", path_length);

            // Print the path:
            print!("Polyline(");
            for (i, p) in path.iter().enumerate() {
                if i > 0 {
                    print!(",");
                }
                print!("({}, {}, {})", p.point.x, p.point.y, p.point.z);
            }
            println!(")");

            // Check that the start and end are correct:
            assert_eq!(path.first().unwrap().point, source_point.point);
            assert_eq!(path.last().unwrap().point, goal_point.point);

            // The path length should be at most as long as the path found by the edge graph.
            assert!(path_length <= edge_graph_distance + 1.0e-6);

            // As a sanity check, it should not be shorter than the straight-line distance.
            assert!(path_length >= distance(&source_point.point, &goal_point.point) - 1.0e-6);
        }
    }
}

// Note: Visualize implementation commented out as it depends on visualization modules
// that are not part of the standalone geometry crate
// impl Visualize for SequenceTree {
//     fn visualize(&self) -> three_d::CpuMesh {
//         line_segments_mesh(self.windows.iter().flat_map(|w| {
//             let origin = w.coverage.origin;
//             let left_exit = w.face_exit.left_point;
//             let right_exit = w.face_exit.right_point;
//             [
//                 (origin.to_three_d(), left_exit.to_three_d()),
//                 (origin.to_three_d(), right_exit.to_three_d()),
//             ]
//         }))
//     }
// }
