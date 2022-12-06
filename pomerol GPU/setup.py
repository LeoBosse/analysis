from setuptools import setup, find_packages
from Cython.Build import cythonize


setup(  name='pomerol',
        version='1.0',
        packages=find_packages(),
        ext_modules=cythonize("cythonized.pyx", annotate=True),
        zip_safe=False)
