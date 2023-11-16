from setuptools import setup, find_packages

setup(
    name='CurveSimilarity',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'scipy'
    ],
    author='Yajing Li/李亚京',
    author_email='liyajing20@mails.tsinghua.edu.cn',
    description='A Python library for calculating the similarity of 1D, 2D, and 3D curves.',
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url='http://github.com/Cynic2019/CurveSimilarity',
    python_requires='>=3.5',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11'
    ]
)
