from setuptools import setup


setup(name='MHcut',
      version='0.1',
      description='Micro-Homology cut finder',
      url='http://github.com/jmonlone/MHcut',
      author='Woltjen & Bourque labs',
      author_email='jean.monlong@gmail.com',
      license='MIT',
      packages=['MHcut'],
      install_requires=['pyfaidx', 'tqdm', 'argparse'],
      entry_points={'console_scripts': ['MHcut = MHcut.run_mhcut:main']})
