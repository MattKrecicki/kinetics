from os.path import join
from glob import glob

import setuptools


DATA_EXTS = {'*.h5'}


def getDataFiles():
    """Return all data files from ``kineitcs/data``"""

    files = []
    for ext in DATA_EXTS:
        files.extend(glob(join('kinetics', 'data', ext)))
    return files


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="kinetics",
    version="0.0.1",
    author="Matt Krecicki",
    author_email="matthewkrecicki@gmail.com",
    description="A reactor kinetics modeling package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    packages=['kinetics', 'kinetics.functions', 'kinetics.data',
              'kinetics.debug', 'kinetics.errors', 'kinetics.examples',
              'kinetics.tests', 'kinetics.containers'],
    package_data={
        'kineitcs.data': ['data/{}'.format(ext) for ext in DATA_EXTS],
    },
    include_package_data=True,
    data_files=[('kinetics/data', getDataFiles()), ],
)
