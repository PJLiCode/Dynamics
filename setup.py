from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        "boltzmann_wealth",  # source file
        ["boltzmann_bindings.cpp"],  # binding file
        include_dirs=[pybind11.get_include()],  # pybind11 directory
        language="c++",  
    ),
]

setup(
    name="boltzmann_wealth", # package name
    ext_modules=ext_modules,
)
