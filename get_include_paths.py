#!/usr/bin/env python3
"""
Get all necessary include paths for C++ development with pybind11
Automatically generates VS Code configuration
"""

import pybind11
import sysconfig
import json
import os
import sys
import platform

def get_include_paths():
    """Get all include paths needed for pybind11 development"""
    paths = [
        # Project includes
        os.path.join(os.getcwd(), "cpp", "include"),
        
        # pybind11 includes
        pybind11.get_include(),
        pybind11.get_include(user=True),
        
        # Python includes
        sysconfig.get_path('include'),
    ]
    
    # Add platform-specific Python includes
    if platform.system() == "Windows":
        # Windows might have includes in Scripts folder too
        python_base = os.path.dirname(sys.executable)
        potential_include = os.path.join(python_base, "..", "include")
        if os.path.exists(potential_include):
            paths.append(os.path.abspath(potential_include))
    
    # Remove duplicates and non-existent paths
    unique_paths = []
    for path in paths:
        normalized_path = os.path.normpath(path)
        if normalized_path not in unique_paths and os.path.exists(normalized_path):
            unique_paths.append(normalized_path)
    
    return unique_paths

def generate_vscode_config():
    """Generate VS Code C++ configuration"""
    paths = get_include_paths()
    
    # Determine the platform-specific configuration
    system = platform.system()
    if system == "Linux":
        config_name = "Linux"
        intellisense_mode = "linux-gcc-x64"
        compiler_path = "/usr/bin/g++"
    elif system == "Darwin":  # macOS
        config_name = "Mac"
        intellisense_mode = "macos-clang-x64"
        compiler_path = "/usr/bin/clang++"
    else:  # Windows
        config_name = "Win32"
        intellisense_mode = "windows-gcc-x64"
        compiler_path = "cl.exe"  # or path to g++ if using MinGW
    
    config = {
        "configurations": [
            {
                "name": config_name,
                "includePath": [
                    "${workspaceFolder}/**"
                ] + paths,
                "defines": [],
                "compilerPath": compiler_path,
                "cStandard": "c11",
                "cppStandard": "c++17",
                "intelliSenseMode": intellisense_mode
            }
        ],
        "version": 4
    }
    
    # Create .vscode directory if it doesn't exist
    os.makedirs(".vscode", exist_ok=True)
    
    # Check if config already exists
    config_path = os.path.join(".vscode", "c_cpp_properties.json")
    if os.path.exists(config_path):
        response = input(f"{config_path} already exists. Overwrite? (y/n): ")
        if response.lower() != 'y':
            print("Configuration not written.")
            return
    
    # Write configuration
    with open(config_path, "w") as f:
        json.dump(config, f, indent=4)
    
    print(f"✓ Created {config_path}")

def test_imports():
    """Test if we can import required modules"""
    print("\nTesting Python module imports:")
    print("-" * 50)
    
    modules = ['pybind11', 'numpy', 'matplotlib', 'streamlit']
    for module in modules:
        try:
            __import__(module)
            print(f"✓ {module:<15} - OK")
        except ImportError:
            print(f"✗ {module:<15} - Not installed")

def main():
    print("Monte Carlo Shielding - Include Path Configuration")
    print("=" * 50)
    print(f"Python executable: {sys.executable}")
    print(f"Python version: {sys.version.split()[0]}")
    print(f"Platform: {platform.system()}")
    print()
    
    print("Include paths for your C++ IDE:")
    print("-" * 50)
    for i, path in enumerate(get_include_paths(), 1):
        print(f"{i}. {path}")
    
    print()
    
    # Generate VS Code config if requested
    if len(sys.argv) > 1 and sys.argv[1] == "--vscode":
        print("Generating VS Code configuration...")
        generate_vscode_config()
    else:
        response = input("Generate VS Code configuration? (y/n): ")
        if response.lower() == 'y':
            generate_vscode_config()
    
    # Test imports
    test_imports()
    
    print("\n" + "=" * 50)
    print("Configuration complete!")
    print("\nNext steps:")
    print("1. Restart VS Code to reload the configuration")
    print("2. Open hello.cpp - the red squiggles should be gone")
    print("3. You should have IntelliSense for pybind11 functions")

if __name__ == "__main__":
    main()