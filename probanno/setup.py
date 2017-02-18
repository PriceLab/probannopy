from setuptools import setup

setup(name='probanno',
      version='0.1',
      description='A stand-alone python package for the probabilistic annotation algorithm',
      author='Terry Farrah<tfarrah@systemsbiology.org>, Brendan King <bking@systemsbiology.org',
      author_email='bking@systemsbiology.org',
      license='MIT',
      packages=['lib', 'test'],
      install_requires=[
          'cobra',
      ],
      zip_safe=False)
