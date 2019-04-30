from setuptools import setup


setup(name='MHcut',
      version='0.1',
      description='Micro-Homology cut finder',
      url='http://github.com/jmonlone/MHcut',
      author='Woltjen & Bourque labs',
      author_email='jean.monlong@gmail.com',
      license='MIT',
      packages=['MHcut'],
      install_requires=['pyfaidx', 'tqdm', 'argparse', 'numpy==1.15.3',
                        'pandas>=0.20.3', 'scikit-learn==0.20.0',
                        'scipy==1.1.0', 'h5py', 'biopython'],
      entry_points={'console_scripts': ['MHcut = MHcut.run_mhcut:main']})
