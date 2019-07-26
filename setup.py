from setuptools import setup, find_packages


setup(name='MHcut',
      version='0.1.1',
      description='Micro-Homology cut finder',
      url='http://github.com/WoltjenLab/MHcut',
      author='Woltjen & Bourque labs',
      author_email='jean.monlong@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=['pyfaidx', 'tqdm', 'argparse', 'numpy==1.15.3',
                        'pandas==0.23.4', 'scikit-learn==0.20.0',
                        'scipy==1.1.0'],
      entry_points={'console_scripts': ['MHcut = MHcut.run_mhcut:main']})
