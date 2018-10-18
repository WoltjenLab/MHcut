from setuptools import setup, find_packages


setup(name='MHcut',
      version='0.1',
      description='Micro-Homology cut finder',
      url='http://github.com/jmonlone/MHcut',
      author='Woltjen & Bourque labs',
      author_email='jean.monlong@gmail.com',
      license='MIT',
      packages=['MHcut'],
      install_requires=find_packages('MHcut'),
      entry_points={'console_scripts': ['MHcut = MHcut.run_mhcut:main']})
