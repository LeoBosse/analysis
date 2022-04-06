from setuptools import setup, find_packages
from Cython.Build import cythonize

setup(  name='pomerol',
        version='1.0',
        packages=find_packages(),
        ext_modules=cythonize([ "main.py",
                                "world.py",
                                "rayleigh_utils.py",
                                "atmosphere.py",
                                "ground_map.py",
                                "multiple_scattering.pyx",
                                "sky_map.py",], annotate=True),
        zip_safe=False,)
