from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(name='MHcut',
      version='1.0.0',
      description='Micro-Homology cut finder',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='http://github.com/WoltjenLab/MHcut',
      author='Woltjen & Bourque labs',
      author_email='jean.monlong@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=['pyfaidx', 'tqdm', 'argparse', 'numpy==1.22.0',
                        'pandas==0.23.4', 'scikit-learn==0.20.0',
                        'scipy==1.1.0'],
      entry_points={'console_scripts': ['MHcut = MHcut.run_mhcut:main']})
