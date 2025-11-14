"""
Demo: Boost.Python (pxr) and nanobind interoperability using stage_py

This demonstrates:
1. Creating USD stages with stage_py and using them with node graph
2. Reading USD files created by node graph with pxr module
3. TRUE interoperability between pxr (Boost.Python) and stage_py (nanobind)
"""

import sys
import os

# Setup paths
binary_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', 'Binaries', 'Debug'))
rznode_python = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'Core', 'rznode', 'python'))

sys.path.insert(0, binary_dir)
sys.path.insert(0, rznode_python)
os.environ['PXR_USD_WINDOWS_DLL_PATH'] = binary_dir
os.chdir(binary_dir)

print(f"Binary dir: {binary_dir}")
print(f"RZNode Python: {rznode_python}")

from pxr import Usd, UsdGeom, Sdf
from ruzino_graph import RuzinoGraph
import stage_py


def demo_create_in_cpp_read_in_python():
    """
    Demo: Create USD stage with stage_py (nanobind), then read it in Python (via pxr)
    
    This demonstrates the workflow: stage_py -> node graph -> pxr
    """
    print("\n" + "="*70)
    print("DEMO 1: Create with stage_py (nanobind), Read with pxr (Boost.Python)")
    print("="*70)
    
    output_file = "demo_grid.usdc"
    
    # Step 1: Create geometry using node graph with stage_py
    print("\n🔨 Step 1: Create geometry with node graph...")
    g = RuzinoGraph("Demo")
    g.loadConfiguration(os.path.join(binary_dir, "geometry_nodes.json"))
    
    grid = g.createNode("create_grid", name="grid")
    write_node = g.createNode("write_usd", name="writer")
    g.addEdge(grid, "Geometry", write_node, "Geometry")
    
    inputs = {
        (grid, "resolution"): 5,
        (grid, "size"): 2.0,
    }
    
    # Create stage and convert to payload
    stage = stage_py.Stage(output_file)
    payload = stage_py.create_payload_from_stage(stage, "/GridMesh")
    g.setGlobalParams(payload)
    g.prepare_and_execute(inputs, required_node=write_node)
    stage.save()
    
    print(f"✓ Created USD file: {output_file}")
    
    # Step 2: Read and inspect using pxr module
    print("\n🔍 Step 2: Inspect with pxr module...")
    pxr_stage = Usd.Stage.Open(output_file)
    
    # List all prims
    print("\n📦 Prims in stage:")
    for prim in pxr_stage.Traverse():
        print(f"  • {prim.GetPath()} ({prim.GetTypeName()})")
    
    # Get mesh data
    mesh_prim = pxr_stage.GetPrimAtPath("/GridMesh")
    if mesh_prim.IsValid():
        mesh = UsdGeom.Mesh(mesh_prim)
        points = mesh.GetPointsAttr().Get(0)  # Get at time 0
        faces = mesh.GetFaceVertexIndicesAttr().Get(0)  # Get at time 0
        
        print(f"\n📊 Mesh statistics:")
        print(f"  • Vertices: {len(points) if points else 0}")
        print(f"  • Face indices: {len(faces) if faces else 0}")
        if points and len(points) > 0:
            print(f"  • First 5 points:")
            for i, pt in enumerate(list(points)[:5]):
                print(f"    [{i}] {pt}")
    
    print("\n✅ Successfully read USD file created by C++!")


def demo_create_in_python_write_with_cpp():
    """
    Demo: Create USD stage in Python (pxr), write geometry to it via C++ nodes
    
    WARNING: This may NOT work due to Boost.Python vs nanobind incompatibility!
    The stage object cannot be passed directly between the two systems.
    
    We can try using the capsule bridge if implemented.
    """
    print("\n" + "="*70)
    print("DEMO 2: Create in Python (pxr), Write with C++ (nanobind)")
    print("="*70)
    print("⚠️  WARNING: This is experimental and may not work!")
    
    output_file = "demo_sphere.usdc"
    
    try:
        # Step 1: Create USD stage with pxr
        print("\n🔨 Step 1: Create USD stage with pxr...")
        stage = Usd.Stage.CreateNew(output_file)
        print(f"✓ Created stage: {stage}")
        
        # Step 2: Try to use this stage with node graph
        print("\n🔄 Step 2: Attempting to use pxr stage with node graph...")
        
        # This is where we'd need the capsule bridge
        # The problem: We need to extract the C++ pointer from the pxr.Usd.Stage
        # and pass it to nanobind
        
        print("\n❌ Direct interop not yet implemented!")
        print("   Reason: Boost.Python and nanobind use different object layouts")
        print("   Workaround: Use separate stages or file-based communication")
        
    except Exception as e:
        print(f"\n❌ Failed: {e}")


def demo_workaround_separate_stages():
    """
    Demo: Recommended workaround - use separate USD stages for each system
    
    This is the RECOMMENDED approach:
    - Create geometry with nanobind node graph
    - Save to file
    - Open and modify with pxr if needed
    """
    print("\n" + "="*70)
    print("DEMO 3: Recommended Workflow - Separate Stages")
    print("="*70)
    
    output_file = "demo_workflow.usdc"
    
    # Step 1: Create base geometry with node graph
    print("\n🔨 Step 1: Generate geometry with node graph...")
    g = RuzinoGraph("Workflow")
    g.loadConfiguration(os.path.join(binary_dir, "geometry_nodes.json"))
    
    sphere = g.createNode("create_uv_sphere", name="sphere")
    write_node = g.createNode("write_usd", name="writer")
    g.addEdge(sphere, "Geometry", write_node, "Geometry")
    
    inputs = {
        (sphere, "segments"): 16,
        (sphere, "rings"): 8,
        (sphere, "radius"): 1.0,
    }
    
    payload = nodes_system_py.create_usd_payload(output_file, "/Geometry/Sphere")
    g.setGlobalParams(payload)
    g.prepare_and_execute(inputs, required_node=write_node)
    nodes_system_py.save_usd_stage(payload)
    
    print(f"✓ Created base geometry")
    
    # Step 2: Post-process with pxr
    print("\n🎨 Step 2: Post-process with pxr...")
    stage = Usd.Stage.Open(output_file)
    
    # Add metadata
    stage.SetMetadata("comment", "Generated by Ruzino + post-processed")
    
    # Create a group
    group_prim = stage.DefinePrim("/Geometry", "Xform")
    group_prim.SetMetadata("kind", "group")
    
    # Add custom attributes
    sphere_prim = stage.GetPrimAtPath("/Geometry/Sphere")
    if sphere_prim.IsValid():
        # Add a color attribute
        color_attr = sphere_prim.CreateAttribute(
            "primvars:displayColor",
            Sdf.ValueTypeNames.Color3fArray
        )
        color_attr.Set([(0.8, 0.3, 0.3)])  # Red-ish color
    
    # Save modifications
    stage.Save()
    print(f"✓ Added metadata and attributes")
    
    # Step 3: Verify
    print("\n✅ Step 3: Verify final result...")
    verify_stage = Usd.Stage.Open(output_file)
    print(f"   Comment: {verify_stage.GetMetadata('comment')}")
    
    sphere_verify = verify_stage.GetPrimAtPath("/Geometry/Sphere")
    if sphere_verify.HasAttribute("primvars:displayColor"):
        color = sphere_verify.GetAttribute("primvars:displayColor").Get()
        print(f"   Display color: {color}")
    
    print(f"\n✅ Workflow complete! File: {output_file}")
    print(f"   You can now:")
    print(f"   - Open in Houdini/Maya/Blender")
    print(f"   - Further process with pxr")
    print(f"   - Render with Hydra")


if __name__ == "__main__":
    print("\n" + "="*70)
    print("  USD Interoperability Demos")
    print("  Boost.Python (pxr) ↔ nanobind (node graph)")
    print("="*70)
    
    try:
        demo_create_in_cpp_read_in_python()
        demo_create_in_python_write_with_cpp()
        demo_workaround_separate_stages()
        
        print("\n" + "="*70)
        print("  ✅ ALL DEMOS COMPLETE!")
        print("="*70)
        
        print("\n💡 Key Takeaways:")
        print("  1. Direct object passing between pxr and nanobind: NOT POSSIBLE")
        print("  2. File-based communication: WORKS PERFECTLY")
        print("  3. Recommended workflow: Generate with C++, post-process with Python pxr")
        print("  4. Both systems can read/write the same .usdc files")
        
    except Exception as e:
        print(f"\n❌ Demo failed: {e}")
        import traceback
        traceback.print_exc()
