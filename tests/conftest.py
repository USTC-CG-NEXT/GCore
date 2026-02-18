"""
pytest configuration for geometry tests
Sets up paths and environment variables
"""

import sys
import os
import platform

# Get binary directory
binary_dir = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "..", "Binaries", "Release"
)
binary_dir = os.path.abspath(binary_dir)

project_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..")
)


sys.path.append(binary_dir)

# os.add_dll_directory is Windows-specific
if platform.system() == 'Windows':
    os.add_dll_directory(project_root + r"\SDK\python")
    os.add_dll_directory(project_root + r"\SDK\OpenUSD\Release\lib")
else:
    # Linux: Add USD Python bindings to path
    usd_python_path = os.path.join(project_root, "SDK", "OpenUSD", "Debug", "lib", "python")
    usd_lib_path = os.path.join(project_root, "SDK", "OpenUSD", "Debug", "lib")
    if os.path.exists(usd_python_path):
        sys.path.insert(0, usd_python_path)
    # Preload USD shared libraries using ctypes (LD_LIBRARY_PATH doesn't work after process starts)
    import ctypes
    import glob
    for so_file in glob.glob(os.path.join(usd_lib_path, "*.so")):
        try:
            ctypes.CDLL(so_file, mode=ctypes.RTLD_GLOBAL)
        except OSError:
            pass  # Some libraries may fail to load, that's OK
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
