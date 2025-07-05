from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="GEddySoft",
    version="4.0.0",
    author="Bernard Heinesch",
    author_email="bernard.heinesch@uliege.be",
    description="A Python package for surface eddy-flux computation, with a specific focus on BVOCs from PTR-TOF-MS",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Bernard-Heinesch/GEddySoft",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
)
