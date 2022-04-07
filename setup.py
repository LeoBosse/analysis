from setuptools import setup, find_packages
from Cython.Build import cythonize

to_cythonize = [        "main.py",
                        "world.py",
                        "rayleigh_utils.py",
                        "atmosphere.py",
                        "ground_map.py",
                        "multiple_scattering.pyx",
                        "sky_map.py"]



setup(  name='pomerol',
        version='1.0',
        packages=find_packages(),
        ext_modules=cythonize(to_cythonize, annotate=True),
        zip_safe=False)
