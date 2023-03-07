#from distutils.core import setup
from setuptools import setup, find_packages
import os
import io

def _read(*parts, **kwargs):
    filepath = os.path.join(os.path.dirname(__file__), *parts)
    encoding = kwargs.pop('encoding', 'utf-8')
    with io.open(filepath, encoding=encoding) as fh:
        text = fh.read()
    return text

def get_requirements(path):
    content = _read(path)
    return [
        req
        for req in content.split("\n")
        if req != '' and not req.startswith('#')
    ]

install_requires = get_requirements('requirements.txt')

setup(
    name='panopticrispr',
    version='0.1',
    author='Samuel C. Markson',
    author_email='smarkson@alum.mit.edu',
    install_requires=install_requires, 
    packages=find_packages(),
    scripts=['scripts/panopticrispr'],
    license='MIT',
    description='for general analysis of single-cell RNA data with CRISPR-Cas9 or similar perturbations',
    include_package_data = True,
    long_description=open('README.rst').read(),
)

