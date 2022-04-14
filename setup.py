import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="flopyrwpt", 
    version="0.0.1",
    author="MARSoluT-ITN",
    author_email="rodrigo.alfonso.perez@upc.edu",
    description="Utilities for configuration of rwpt models with flopy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/models7/flopyrwpt",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License", 
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'flopy >= 3.3.5',
    ]
)
