
from setuptools import setup, find_packages

setup(
    author="Shyam Katnagallu",
    author_email='s.katnagallu@mpie.de',
    python_requires='>=3.8',
    
    description="Advanced analysi of iron oxides with correlation histograms",
    install_requires = ['numpy', 'matplotlib', 'pandas', 'scipy'],
            
    include_package_data=True,
    keywords='CorrelationHistograms, AtomProbe',
    name='CorrelationHistogram',
    packages=find_packages(include=['correlation_histogram', 'dissociation_tracks.*']),
)