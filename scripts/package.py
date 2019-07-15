#Copyright 2014-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install vizaly-foresight
#
# You can edit this file again by typing:
#
#     spack edit vizaly-foresight
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *
import os

class VizalyForesight(Package):
    """Foresight is a tool that allows users to check how different compression algorithsm and compression parameters will affect their sumulation data."""

    homepage = "https://github.com/lanl/VizAly-Foresight"
    url      = "https://github.com/lanl/VizAly-Foresight/archive/master.tar.gz"
    git      = "https://github.com/lanl/VizAly-Foresight.git"

    #version('ldav', branch='master')
    version('latest',branch='cmake_install')

    # CBench Depends
    depends_on('mpi')
    depends_on('c-blosc')
    depends_on('sz')
    depends_on('zfp')
    depends_on('hdf5')

    # PAT Depends
    #depends_on('python')
    #depends_on('py-pandas')
    #depends_on('py-numpy')
    #depends_on('py-argparse')
    #depends_on('libpng')

    # Build System
    depends_on('cmake@3.10.2', type='build')

    root_cmakelists_dir = 'CBench'
    #build_directory = 'build-linux' #this does nothing

    def cmake_args(self):
        spec = self.spec
        #working_dir('build-linux', create=True) #this does nothing

        args = []
        args.append('-DWORKING_DIRECTORY=CBench')
        #args.append('-DCMAKE_BUILD_TYPE=DEBUG')

        if 'c-blosc' in spec:
            args.append('-DCBENCH_ENABLE_BLOSC=ON')
            #cblosc_libs = spec['c-blosc'].libs.joined(';')
            shared = True if '+shared' in self.spec else False
            cblosc_libs = find_libraries('libblosc*', root=spec['c-blosc'].prefix, shared=shared, recursive=True).joined(';')
            lz4_libs = spec['lz4'].libs.joined(';')
            snappy_libs = spec['snappy'].libs.joined(';')
            args.append('-DBLOSC_LIBRARY={0}'.format(cblosc_libs+';'+lz4_libs+';'+snappy_libs))
            cblosc_incl = spec['c-blosc'].prefix.include
            args.append('-DBLOSC_INCLUDE_PATH={0}'.format(cblosc_incl))
        if 'zfp' in spec:
            args.append('-DCBENCH_ENABLE_ZFP=ON')
            zfp_libs = spec['zfp'].libs.joined(';')
            args.append('-DZFP_LIBRARY={0}'.format(zfp_libs))
            zfp_incl = spec['zfp'].prefix.include
            args.append('-DZFP_INCLUDE_PATH={0}'.format(zfp_incl))
        if 'sz' in spec:
            args.append('-DCBENCH_ENABLE_SZ=ON')
            shared = True if '+shared' in self.spec else False
            #sz_libs = spec['sz'].libs.joined(';')
            sz_libs = find_libraries('libSZ*', root=spec['sz'].prefix, shared=shared, recursive=True).joined(';')
            args.append('-DSZ_LIBRARY={0}'.format(sz_libs))
            sz_incl = spec['sz'].prefix.include
            args.append('-DSZ_INCLUDE_PATH={0}'.format(sz_incl))

            z_libs = find_libraries('libzlib*', root=spec['sz'].prefix, shared=shared, recursive=True).joined(';')
            args.append('-DZLIB_LIBRARY={0}'.format(z_libs))
            zstd_libs = find_libraries('libzstd*', root=spec['sz'].prefix, shared=shared, recursive=True).joined(';')
            args.append('-DZSTD_LIBRARY={0}'.format(zstd_libs))
        if 'hdf5' in spec:
            args.append('-DCBENCH_ENABLE_NYX_LOADER=ON')
            hdf5_libs = spec['hdf5'].libs.joined(';')
            args.append('-DHDF5_LIBRARY={0}'.format(hdf5_libs))
            hdf5_incl = spec['hdf5'].prefix.include
            args.append('-DHDF5_INCLUDE_PATH={0}'.format(hdf5_incl))

        return args

    def install(self, spec, prefix):
        os.rename('CBench', 'src')
        cmake('src',*std_cmake_args + self.cmake_args())
        make('install')
