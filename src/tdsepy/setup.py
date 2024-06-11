from setuptools import setup

setup(
    name='tdsepy',
    version='0.0.0',
    description='A package for solving the Time-Dependent Schrödinger Equation (TDSE)',
    author='Joshua Mann',
    author_email='jomann@ucla.edu',
    packages=["tdsepy", "tdsepy.helpers"],
    package_data={"tdsepy": ["tdsepy/*.pyi", "py.typed"]},
    install_requires=[
        'numpy',
        'scipy',
    ],
)