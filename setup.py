from setuptools import find_packages, setup


setup(
    name="mt_ckd",
    version="3.5.0",
    author="R. Menzel, R. Pernak",
    author_email="",
    description="MT-CKD continuum model",
    url="",
    python_requires=">=3.5",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: ",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "matplotlib",
        "netCDF4",
        "numpy",
    ],
)
