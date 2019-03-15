# use CentOS 7 as base image
FROM centos:7

# install packages
RUN yum install -y epel-release && \
    yum install -y doxygen git gsl-devel hdf5-devel help2man make \
                   perl-CPAN perl-Thread-Queue perl-devel \
                   openssl-devel selinux-policy sudo wget

## install gcc
RUN yum install -y centos-release-scl && \
    yum install -y devtoolset-6-gcc devtoolset-6-gcc-c++ \
                   devtoolset-6-libquadmath-devel devtoolset-6-gcc-gfortran

# install OpenCL
RUN yum install -y ocl-icd-devel

# install CMake
WORKDIR /src
RUN source /opt/rh/devtoolset-6/enable && \
    wget https://cmake.org/files/v3.11/cmake-3.11.0.tar.gz && \
    tar -zxvf cmake-3.11.0.tar.gz && cd cmake-3.11.0 && \
    ./bootstrap --prefix=/usr/local && \
    make -j $(getconf _NPROCESSORS_ONLN) install

# install OpenMPI
ENV MPI_URL https://www.open-mpi.org/software/ompi/v2.1/downloads
ENV MPI_VERSION 2.1.3
RUN wget -nv ${MPI_URL}/openmpi-${MPI_VERSION}.tar.gz
RUN tar -xf openmpi-${MPI_VERSION}.tar.gz
RUN source /opt/rh/devtoolset-6/enable && cd openmpi-${MPI_VERSION} && \
    CFLAGS=-O3 && \
    CXXFLAGS=-O3 && \
    ./configure --prefix=/usr \
                --sysconfdir=/mnt/0 \
                --with-slurm=no && \
    make -j $(getconf _NPROCESSORS_ONLN) install
RUN echo 'plm_rsh_agent = false' >> /mnt/0/openmpi-mca-params.conf

# install M4
RUN source /opt/rh/devtoolset-6/enable && \
    wget -O m4-1.4.9.tar.gz https://mirrors.ocf.berkeley.edu/gnu/m4/m4-1.4.9.tar.gz && \
    tar -zvxf m4-1.4.9.tar.gz && \
    cd m4-1.4.9 && \
    ./configure --prefix=/usr && \
    make -j $(getconf _NPROCESSORS_ONLN) install

# install Autoconf
RUN source /opt/rh/devtoolset-6/enable && \
    wget https://mirrors.ocf.berkeley.edu/gnu/autoconf/autoconf-2.69.tar.gz && \
    gunzip autoconf-2.69.tar.gz && \
    tar -xvf autoconf-2.69.tar && \
    cd autoconf-2.69 && \
    ./configure --prefix=/usr && \
    make -j $(getconf _NPROCESSORS_ONLN) install

# install AutoMake
RUN source /opt/rh/devtoolset-6/enable && \
    wget https://mirrors.ocf.berkeley.edu/gnu/automake/automake-1.15.tar.gz && \
    tar -xvzf automake-1.15.tar.gz && \
    cd automake-1.15 && \
    ./configure --prefix=/usr && \
    make -j $(getconf _NPROCESSORS_ONLN) install

# create env script
RUN echo "source /opt/rh/devtoolset-6/enable" > env.sh && \
    echo "ldconfig" >> env.sh && \
    echo "export OpenCL_LIBRARY=/usr/lib64/libOpenCL.so" >> env.sh
