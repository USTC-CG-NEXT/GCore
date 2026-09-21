"""
pytest configuration for geometry tests
Sets up paths and environment variables
"""

import sys
import os

project_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..")
)

# RZ_BUILD_TYPE lets this suite run against Binaries/Debug (native debugging
# with PDBs); fall back to Release when that build type isn't built. Mirrors
# source/tests/conftest.py.
_build_type = os.environ.get("RZ_BUILD_TYPE", "Release")
binary_dir = os.path.abspath(os.path.join(project_root, "Binaries", _build_type))
if not os.path.isdir(binary_dir):
    binary_dir = os.path.abspath(os.path.join(project_root, "Binaries", "Release"))

# All test outputs (.usdc) must live under Binaries/Release/test_output/geometry/.
# Binaries/ is gitignored at the repo root, so nothing here enters the source
# tree. Tests import this via `from conftest import OUTPUT_DIR`.
OUTPUT_DIR = os.path.join(binary_dir, "test_output", "geometry")
os.makedirs(OUTPUT_DIR, exist_ok=True)

sys.path.append(binary_dir)
if sys.platform == "win32":
    os.add_dll_directory(project_root + r"\SDK\python")
    os.add_dll_directory(project_root + r"\SDK\OpenUSD\Release\lib")
# Set PXR_USD_WINDOWS_DLL_PATH so USD can find its DLLs
os.environ["PXR_USD_WINDOWS_DLL_PATH"] = binary_dir
print(f"Set PXR_USD_WINDOWS_DLL_PATH={binary_dir}")

# Add to Python path
sys.path.insert(0, binary_dir)

# Add rznode python path
rznode_python = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "Core", "rznode", "python"
)
sys.path.insert(0, os.path.abspath(rznode_python))

# Change to binary dir so DLLs can be loaded
os.chdir(binary_dir)
print(f"Changed working directory to: {os.getcwd()}")
