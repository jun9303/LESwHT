import os
import math
import datetime as dt
import numpy as np
import trimesh

class DimplePatternGenerator:
    def __init__(self, baseline_dir):
        self.baseline_dir = baseline_dir
        
        self.topology_map = {
            0: ('case1_circular.stl', False),
            1: ('case2_ngCricular.stl', False),
            2: ('case3_ngTearDrop.stl', False),
            3: ('case3_ngTearDrop.stl', True),  
            4: ('case4_ngDiamond.stl', False),
            5: ('case5_isosceles.stl', False),
            6: ('case5_isosceles.stl', True)    
        }
        
        self.topology_descriptions = {
            0: "Cylinder: Circular top with cylindrical carve.",
            1: "Spherical: Circular top with smooth depth variation.",
            2: "Teardrop Down: Teardrop shape (isosceles + circle at foot).",
            3: "Teardrop Up: X-flipped Teardrop Down.",
            4: "Diamond: Diamond-like shape (overlapping teardrops).",
            5: "Triangle Down: Isosceles shape with linear varying depth.",
            6: "Triangle Up: X-flipped Triangle Down."
        }

    def _streamwise_spacing(self, topo_idx):
        """Return base streamwise center-to-center spacing for each topology family."""
        if topo_idx in [0, 1]:
            return 5.0 * (math.sqrt(3) / 2.0)
        if topo_idx in [2, 3, 4]:
            return 2.5 * (4.0 - math.sqrt(3))
        if topo_idx in [5, 6]:
            return 5.0
        return 5.0

    def _cleanup_mesh(self, mesh):
        """Apply compatibility-safe cleanup for different trimesh versions."""
        if hasattr(mesh, "remove_duplicate_faces"):
            mesh.remove_duplicate_faces()

        mesh.remove_unreferenced_vertices()
        mesh.fix_normals()

    def _boolean_difference(self, block, all_punches):
        """Robust boolean subtraction using manifold engine when available."""
        try:
            return trimesh.boolean.difference(
                [block, all_punches],
                engine='manifold',
                check_volume=False
            )
        except Exception:
            # Fallback to trimesh default behavior if manifold backend is unavailable.
            return block.difference(all_punches)

    def _align_streamwise_wallnormal_spanwise(self, mesh):
        """Keep coordinates in (x, y, z) = (streamwise, wall-normal, spanwise), with top at y=0."""
        ymax = float(np.max(mesh.vertices[:, 1]))
        if abs(ymax) > 1.0e-12:
            mesh.vertices[:, 1] -= ymax

    def _effective_scale_for_fill(self, scale_val):
        """Avoid exact tangency for fill booleans when user requests scale=1.0."""
        if np.isclose(scale_val, 1.0):
            tol = 1.0e-3
            return 1.0 - tol
        return scale_val

    def _create_cylindrical_baseline(self, target_depth, scale_val, segments=128):
        """Construct a tapered circular dimple surface with top boundary at y=0."""
        top_radius = 2.45 * scale_val
        bottom_radius = 0.95 * top_radius
        angles = np.linspace(0.0, 2.0 * math.pi, segments, endpoint=False)

        x_top = top_radius * np.cos(angles)
        z_top = top_radius * np.sin(angles)
        x_bottom = bottom_radius * np.cos(angles)
        z_bottom = bottom_radius * np.sin(angles)

        top_ring = np.column_stack([x_top, np.zeros(segments), z_top])
        bottom_ring = np.column_stack([x_bottom, np.full(segments, target_depth), z_bottom])
        bottom_center = np.array([[0.0, target_depth, 0.0]])

        vertices = np.vstack([top_ring, bottom_ring, bottom_center])
        center_idx = 2 * segments

        faces = []
        for i in range(segments):
            j = (i + 1) % segments

            top_i = i
            top_j = j
            bottom_i = segments + i
            bottom_j = segments + j

            # Side wall.
            faces.append([top_i, bottom_i, bottom_j])
            faces.append([top_i, bottom_j, top_j])

            # Bottom disk.
            faces.append([bottom_i, center_idx, bottom_j])

        baseline = trimesh.Trimesh(vertices=vertices, faces=np.array(faces), process=False)
        baseline.fix_normals()
        return baseline

    def _top_boundary_xz(self, mesh):
        """Return 2D (x, z) points of the open top boundary at maximum y."""
        edges = mesh.edges
        edges_sorted = np.sort(edges, axis=1)
        _, inverse, counts = np.unique(edges_sorted, axis=0, return_inverse=True, return_counts=True)
        boundary_edge_indices = np.where(counts[inverse] == 1)[0]
        boundary_edges = edges[boundary_edge_indices]
        boundary_vertices = np.unique(boundary_edges.reshape(-1))
        vb = mesh.vertices[boundary_vertices]

        y_top = np.max(vb[:, 1])
        top = vb[np.abs(vb[:, 1] - y_top) < 1.0e-8][:, [0, 2]]
        return top

    def _anchor_teardrop(self, mesh):
        """Anchor teardrop at the center of its (half-)circular region."""
        top = self._top_boundary_xz(mesh)
        if len(top) == 0:
            return 0.0, 0.0

        z_min = np.min(top[:, 1])
        z_max = np.max(top[:, 1])
        z_span = max(z_max - z_min, 1.0e-12)
        z_tol = 0.01 * z_span

        low = top[top[:, 1] <= z_min + z_tol]
        high = top[top[:, 1] >= z_max - z_tol]

        if len(low) == 0 or len(high) == 0:
            # Fallback to box-center if numerical filtering misses extrema points.
            return 0.5 * (np.min(top[:, 0]) + np.max(top[:, 0])), 0.5 * (z_min + z_max)

        x_low = np.mean(low[:, 0])
        x_high = np.mean(high[:, 0])
        return 0.5 * (x_low + x_high), 0.5 * (z_min + z_max)

    def _anchor_isosceles_foot_midpoint(self, mesh):
        """Anchor isosceles at midpoint of its foot edge."""
        top = self._top_boundary_xz(mesh)
        if len(top) == 0:
            return 0.0, 0.0

        x_min = np.min(top[:, 0])
        x_max = np.max(top[:, 0])
        x_span = max(x_max - x_min, 1.0e-12)
        x_tol = 0.02 * x_span

        left = top[top[:, 0] <= x_min + x_tol]
        right = top[top[:, 0] >= x_max - x_tol]

        def span_z(points):
            if len(points) == 0:
                return -1.0
            return float(np.max(points[:, 1]) - np.min(points[:, 1]))

        use = left if span_z(left) >= span_z(right) else right
        if len(use) == 0:
            return 0.5 * (x_min + x_max), 0.0

        x_anchor = np.mean(use[:, 0])
        z_anchor = 0.5 * (np.min(use[:, 1]) + np.max(use[:, 1]))
        return x_anchor, z_anchor

    def process_baseline_into_tool(self, topo_idx, depth_val, scale_val, stream_sx=1.0):
        """Loads/scales the dimple and converts it into a punch solid for booleans."""
        target_depth = -(depth_val / 180.0)

        if topo_idx == 0:
            # Case 0 baseline STL has disconnected/non-manifold pieces; generate it analytically.
            # Build unscaled in-plane geometry, then apply common scaling path below.
            base_mesh = self._create_cylindrical_baseline(target_depth, 1.0)
        else:
            filename, is_flipped = self.topology_map[topo_idx]
            filepath = os.path.join(self.baseline_dir, filename)

            # Load mesh and merge any duplicate vertices
            base_mesh = trimesh.load_mesh(filepath)
            base_mesh.process()

            # 1. Flip X if required
            if is_flipped:
                base_mesh.vertices[:, 0] *= -1.0
                base_mesh.invert() # Fix normal winding after mirror

            # 2. Translate so the top boundary is exactly at y=0
            max_y = np.max(base_mesh.vertices[:, 1])
            base_mesh.vertices[:, 1] -= max_y

            # 3. Depth Scaling
            current_depth = np.max(base_mesh.vertices[:, 1]) - np.min(base_mesh.vertices[:, 1])
            d_scale = abs(target_depth / current_depth) if current_depth != 0 else 1.0
            base_mesh.vertices[:, 1] *= d_scale

            # 5. Topology-specific in-plane anchor alignment before tiling
            if topo_idx in [2, 3]:
                ax, az = self._anchor_teardrop(base_mesh)
                base_mesh.vertices[:, 0] -= ax
                base_mesh.vertices[:, 2] -= az
            elif topo_idx in [5, 6]:
                ax, az = self._anchor_isosceles_foot_midpoint(base_mesh)
                base_mesh.vertices[:, 0] -= ax
                base_mesh.vertices[:, 2] -= az

        # Apply streamwise stretch (x only), then isotropic pattern scaling (x and z).
        # Final in-plane multipliers: x -> stream_sx * scale_val, z -> scale_val.
        base_mesh.vertices[:, 0] *= stream_sx
        base_mesh.vertices[:, 0] *= scale_val
        base_mesh.vertices[:, 2] *= scale_val

        # Convert open dimple surface into a subtraction punch by capping boundary to a top apex.
        edges = base_mesh.edges
        edges_sorted = np.sort(edges, axis=1)
        _, inverse, counts = np.unique(edges_sorted, axis=0, return_inverse=True, return_counts=True)
        boundary_edge_indices = np.where(counts[inverse] == 1)[0]
        boundary_edges = edges[boundary_edge_indices]

        top_vertex = np.array([[0.0, 5.0, 0.0]])
        new_vertices = np.vstack([base_mesh.vertices, top_vertex])
        top_v_idx = len(new_vertices) - 1

        new_faces = list(base_mesh.faces)
        for edge in boundary_edges:
            new_faces.append([edge[1], edge[0], top_v_idx])

        punch = trimesh.Trimesh(vertices=new_vertices, faces=new_faces, process=False)
        self._cleanup_mesh(punch)
        return punch, abs(target_depth)
    
    def get_centers(self, topo_idx, sx_scale, compact=False):
        """Calculates the (X, Z) coordinates for the 3-repetition pattern."""
        centers = []

        dx = self._streamwise_spacing(topo_idx) * sx_scale

        if not compact:
            grid_range = range(-3, 4)
            for r in grid_range:
                for c in grid_range:
                    z = r * 5.0
                    if c % 2 != 0:
                        z += 2.5
                    x = c * dx
                    centers.append((x, z))
            return centers

        # Compact periodic cell:
        # - fixed span in z: -5 <= z <= 5
        # - x span contains 2 upstream + 2 downstream spacings from the center
        # Ghost punches are included around the compact bounds to carve periodic edges cleanly.
        x_min = -2.0 * dx
        x_max = 2.0 * dx
        z_min = -5.0
        z_max = 5.0
        ghost = 2.6

        for r in range(-4, 5):
            for c in range(-8, 9):
                z = r * 5.0
                if c % 2 != 0:
                    z += 2.5
                x = c * dx
                if (x_min - ghost) <= x <= (x_max + ghost) and (z_min - ghost) <= z <= (z_max + ghost):
                    centers.append((x, z))
                
        return centers

    def generate(self, topo_idx, depth_val, scale_val, domain_sx, compact=False):
        print("="*60)
        print("GENERATING DIMPLED SLAB MESH:")
        print(f"Topology ID    : {topo_idx}")
        print(f"Target Depth   : {depth_val} wall units")
        if np.isclose(scale_val, 1.0):
            print(f"Input Scale    : {scale_val}")
        print("="*60)
        
        # 1. Get the transformed baseline dimple punch and pattern centers
        eff_scale = self._effective_scale_for_fill(scale_val)
        dimple_punch, actual_depth = self.process_baseline_into_tool(
            topo_idx,
            depth_val,
            eff_scale,
            stream_sx=domain_sx
        )
        centers = self.get_centers(topo_idx, domain_sx, compact=compact)
        
        # 2. Place all dimple surfaces
        print("Arranging dimple pattern...")
        dimples = []
        for (cx, cz) in centers:
            d = dimple_punch.copy()
            d.apply_translation([cx, 0.0, cz])
            dimples.append(d)

        # Combine dimples into one subtraction tool
        all_dimples = trimesh.util.concatenate(dimples)

        # 3. Create base block whose top is at y=0 for non-dimple flat fill
        print("Constructing solid bounding volume...")
        if compact:
            dx = self._streamwise_spacing(topo_idx) * domain_sx
            min_x = -2.0 * dx
            max_x = 2.0 * dx
            min_z = -5.0
            max_z = 5.0
            print(f"Compact periodic domain: x in [{min_x:.6f}, {max_x:.6f}], z in [{min_z:.6f}, {max_z:.6f}]")
        else:
            min_x = min([c[0] for c in centers]) - 5.0
            max_x = max([c[0] for c in centers]) + 5.0
            min_z = min([c[1] for c in centers]) - 5.0
            max_z = max([c[1] for c in centers]) + 5.0

        width = max_x - min_x
        length = max_z - min_z
        block_depth = actual_depth + 2.0

        block = trimesh.creation.box(extents=[width, block_depth, length])
        block.apply_translation([
            (max_x + min_x) / 2.0,
            -block_depth / 2.0,
            (max_z + min_z) / 2.0
        ])

        # 4. Carve dimples from the block so top-plane cavities are filled at y=0.
        print("Performing mathematical boolean difference (Carving holes)...")
        final_mesh = self._boolean_difference(block, all_dimples)
        if final_mesh is None:
            raise RuntimeError("Boolean difference failed to return a mesh.")

        self._cleanup_mesh(final_mesh)
        self._align_streamwise_wallnormal_spanwise(final_mesh)
        output_mesh = final_mesh

        # 6. Export
        filename = f"{topo_idx}_{depth_val}_{scale_val}_{domain_sx}.stl"
        output_mesh.export(filename)
        
        print(f"Success! STL saved as: {filename}")
        print(f"Slab watertight: {output_mesh.is_watertight}")
        print(f"Mesh components: {len(output_mesh.split(only_watertight=False))}\n")

    def generate_to_path(self, topo_idx, depth_val, scale_val, domain_sx, output_path, compact=True):
        """Generate slab and export as ASCII STL to a requested path."""
        print("=" * 60)
        print("GENERATING DIMPLED SLAB MESH:")
        print(f"Topology ID    : {topo_idx}")
        print(f"Target Depth   : {depth_val} wall units")
        if np.isclose(scale_val, 1.0):
            print(f"Input Scale    : {scale_val}")
        print("=" * 60)

        eff_scale = self._effective_scale_for_fill(scale_val)
        dimple_punch, actual_depth = self.process_baseline_into_tool(
            topo_idx,
            depth_val,
            eff_scale,
            stream_sx=domain_sx
        )
        centers = self.get_centers(topo_idx, domain_sx, compact=compact)

        print("Arranging dimple pattern...")
        dimples = []
        for (cx, cz) in centers:
            d = dimple_punch.copy()
            d.apply_translation([cx, 0.0, cz])
            dimples.append(d)

        all_dimples = trimesh.util.concatenate(dimples)

        print("Constructing solid bounding volume...")
        if compact:
            dx = self._streamwise_spacing(topo_idx) * domain_sx
            min_x = -2.0 * dx
            max_x = 2.0 * dx
            min_z = -5.0
            max_z = 5.0
            print(f"Compact periodic domain: x in [{min_x:.6f}, {max_x:.6f}], z in [{min_z:.6f}, {max_z:.6f}]")
        else:
            min_x = min([c[0] for c in centers]) - 5.0
            max_x = max([c[0] for c in centers]) + 5.0
            min_z = min([c[1] for c in centers]) - 5.0
            max_z = max([c[1] for c in centers]) + 5.0

        width = max_x - min_x
        length = max_z - min_z
        block_depth = actual_depth + 2.0

        block = trimesh.creation.box(extents=[width, block_depth, length])
        block.apply_translation([
            (max_x + min_x) / 2.0,
            -block_depth / 2.0,
            (max_z + min_z) / 2.0
        ])

        print("Performing mathematical boolean difference (Carving holes)...")
        final_mesh = self._boolean_difference(block, all_dimples)
        if final_mesh is None:
            raise RuntimeError("Boolean difference failed to return a mesh.")

        self._cleanup_mesh(final_mesh)
        self._align_streamwise_wallnormal_spanwise(final_mesh)

        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        final_mesh.export(output_path, file_type='stl_ascii')

        print(f"Success! STL saved as: {output_path}")
        print(f"Slab watertight: {final_mesh.is_watertight}")
        print(f"Mesh components: {len(final_mesh.split(only_watertight=False))}\n")


def case_name_from_params(topology, depth, scale, stretch, when=None):
    if when is None:
        when = dt.datetime.now()
    stamp = when.strftime("%Y%m%d_%H%M%S")
    return f"{topology}_{depth}_{scale}_{stretch}_{stamp}"


def prepare_case_directories(case_root):
    os.makedirs(case_root, exist_ok=True)
    subdirs = [
        "ftr",
        "field",
        "field_avg",
        "post_inst",
        "post_avg",
        "grid",
        "ibmpre",
        "geometry",
    ]
    for subdir in subdirs:
        os.makedirs(os.path.join(case_root, subdir), exist_ok=True)


def write_case_env(case_root, case_name, stl_path):
    env_path = os.path.join(os.path.dirname(case_root), "running_case.env")
    with open(env_path, "w", encoding="utf-8") as f:
        f.write("# Active LESwHT run context. You may edit this file to switch active case.\n")
        f.write(f"LESWHT_CASE_NAME={case_name}\n")
        f.write(f"LESWHT_OUTPUT_ROOT={case_root}\n")
        f.write(f"LESWHT_GEOMETRY_STL={stl_path}\n")
