from setuptools import setup, find_packages

setup(
    name='NetDes',
    version='1.0.1',
    packages=['NetDes'],  # Automatically find all packages in the directory
    description='Network modeling based on the Dynamic equation simulation (NetDes)',
    license='MIT',
    install_requires=[
        'numpy>=1.26.4',         # Specify the version for numpy
        'torch>=1.12.1',         # Specify the version for PyTorch
        'matplotlib>=3.5.1',     # Specify the version for matplotlib
        'pandas>=2.2.2',         # Specify the version for pandas
        'scipy>=1.13.1',         # Specify the version for scipy
        'scikit-learn>=1.0.2'    # Specify the version for scikit-learn
    ],
    author='Yukai You',
    author_email='you.yu@northeastern.edu',
    url='https://github.com/lusystemsbio/NetDes',
    download_url='https://github.com/lusystemsbio/NetDes/archive/refs/heads/main.zip',
    keywords=['NetDes', 'TF', 'ODE', 'GRN', 'trajectories'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.9',  # Ensure it works for Python 3.9 and above
)
