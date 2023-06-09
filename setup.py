# Copyright (C) 2021-2023 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

from setuptools import setup

setup(
    name='kegg2bipartitegraph',
    url='https://github.com/AuReMe/esmecata',
    license='GPLv3+',
    description=
    'Reconstruct metabolic bipartite graph using KEGG',
    author='AuReMe',
    author_email='gem-aureme@inria.fr',
    packages=['kegg2bipartitegraph'],
    package_dir={'kegg2bipartitegraph': 'kegg2bipartitegraph'},
    entry_points={
        'console_scripts': [
            'kegg2bipartitegraph = kegg2bipartitegraph.__main__:main',
        ]
    },
    install_requires=['bioservices', 'python-libsbml', 'requests', 'networkx'],
)