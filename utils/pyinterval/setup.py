#! /usr/bin/env python

# Copyright (c) 2008-2016 Stefano Taschini <taschini@ieee.org>
# All rights reserved.
# See LICENSE for details.

from setuptools import setup


def read_long_description():
    import io
    parts = []
    for filename in 'README.rst', 'CHANGES.rst':
        with io.open(filename, encoding='utf-8') as f:
            parts.append(f.read())
    return '\n'.join(parts)


metadata = dict(
    author           = 'Stefano Taschini',
    author_email     = 'taschini@gmail.com',
    description      = 'Interval arithmetic in Python',
    install_requires = ['crlibm>=1.0.3,==1.*', 'six>=1.10'],
    keywords         = 'interval crlibm',
    license          = 'BSD',
    long_description = read_long_description(),
    name             = 'pyinterval',
    package_data     = dict(interval=['LICENSE']),
    packages         = ['interval'],
    platforms        = '',
    test_suite       = 'test',
    tests_require    = [],
    url              = "https://github.com/taschini/pyinterval",
    version          = '1.2.1.dev0',
    zip_safe         = False,
    extras_require   = {
        'develop': [
            'Sphinx',
            'flake8',
            'pytest-cov',
            'tox',
            'zest.releaser[recommended]',
        ]
    },
    classifiers      = [
        # A subset of http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Mathematics'
    ],
)


def check_pip():
    try:
        from pip._vendor.pkg_resources import parse_requirements
    except ImportError:
        return

    reqs = metadata['install_requires']

    try:
        [r for s in reqs for r in parse_requirements(s)]
    except ValueError:
        exit("\n"
             "To install this library your must upgrade pip.\n"
             "Please run:\n"
             "   python -m pip install --upgrade pip\n")


if __name__ == '__main__':
    check_pip()
    setup(**metadata)
