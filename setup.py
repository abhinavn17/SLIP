from setuptools import setup, find_packages

# Read the requirements from the requirements.txt file
with open('requirements.txt') as f:
    requirements = f.read().splitlines()


setup(
    name='SLIP',  
    version='2.1.0',  
    author='Abhinav Narayan',
    description='Spectral Line Imaging Pipeline for GMRT and uGMRT data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    packages=find_packages(),    
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'slip = slip.__main__:main',
            'slip_config = slip.__main__:make_config_file',
            'cube_process = slip.cube_process:main',
            'chan_finder = slip.chan_finder:main',
        ]
    },
    include_package_data=True
)
