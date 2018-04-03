from setuptools import setup


with open('oospectro/__init__.py') as fid:
    for line in fid:
        if line.startswith('__version__'):
            VERSION = line.strip().split()[-1][1:-1]
            break


setup(
    name='oospectro',
    version=VERSION,
    url='https://github.com/sciunto-org/oospectro',
    maintainer='F. Boulogne',
    maintainer_email='devel@sciunto.org',
    license='BSD',
    long_description=open('README.md').read(),
    include_package_data=True,
    platforms=['all'],
    description='A library that calculates thicknesses from ocean optics spectra.',
    packages=['oospectro',],
    provides=['oospectro',],
    install_requires = ['numpy', 'scipy', 'matplotlib', 'scikit-image']
    )
