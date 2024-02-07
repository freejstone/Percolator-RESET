import setuptools
import versioneer

with open("README.md", "r") as fh:
    long_description = fh.read()

exec(open('percolator_RESET/_version.py').read())

setuptools.setup(
    name="Percolator RESET",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Jack Freestone",
    author_email="jfre0619@uni.sydney.edu.au",
    description="Percolator RESET (REScoring via Estimating and Training)",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/freejstone/Percolator-RESET",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    python_requires='>=3.9, <=3.11',
    install_requires=[
        'pandas',
        'numpy',
        'scikit-learn>=1.4'
    ],
)