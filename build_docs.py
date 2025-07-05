#!/usr/bin/env python
"""
Build the documentation for GEddySoft.
Run this script from the project root directory.
"""

import os
import sys
import subprocess
import webbrowser
from pathlib import Path

def build_docs():
    """Build the Sphinx documentation."""
    docs_dir = Path(__file__).parent / 'docs'
    build_dir = docs_dir / 'build' / 'html'
    
    # Clean build directory
    if build_dir.exists():
        import shutil
        shutil.rmtree(build_dir)
    
    # Ensure docs/build directory exists
    build_dir.parent.mkdir(parents=True, exist_ok=True)
    
    # Build documentation
    if sys.platform.startswith('win'):
        # Windows
        result = subprocess.run(['cmd', '/c', 'make.bat', 'html'], 
                             cwd=docs_dir, 
                             capture_output=True, 
                             text=True)
    else:
        # Unix-like systems
        result = subprocess.run(['make', 'html'], 
                             cwd=docs_dir, 
                             capture_output=True, 
                             text=True)
    
    if result.returncode != 0:
        print("Error building documentation:")
        print(result.stderr)
        return False
    
    # Open in browser if build successful
    index_path = build_dir / 'index.html'
    if index_path.exists():
        print(f"\nDocumentation built successfully!")
        print(f"Opening {index_path} in your browser...")
        webbrowser.open(index_path.as_uri())
    else:
        print("Documentation build failed: no index.html found")
        return False
    
    return True

if __name__ == '__main__':
    build_docs()
