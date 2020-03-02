# use CentOS 7 as base image
FROM centos:7

# install packages
RUN yum install -y epel-release && \
    yum install -y doxygen git gsl-devel hdf5 hdf5-static \
                   hdf5-openmpi-devel hdf5-openmpi-static \
                   help2man make perl-CPAN perl-Thread-Queue \
                   perl-devel python36-devel python36-pip openmpi \
                   openmpi-devel openssl-devel selinux-policy \
                   sqlite-devel sudo wget

# install gcc
RUN yum install -y centos-release-scl && \
    yum install -y devtoolset-6-gcc devtoolset-6-gcc-c++ \
                   devtoolset-6-libquadmath-devel devtoolset-6-gcc-gfortran

# install OpenCL
RUN yum install -y ocl-icd-devel

# install CMake
WORKDIR /src
RUN source /opt/rh/devtoolset-6/enable && \
    wget https://cmake.org/files/v3.14/cmake-3.14.0.tar.gz && \
    tar -zxvf cmake-3.14.0.tar.gz && cd cmake-3.14.0 && \
    ./bootstrap --prefix=/usr/local && \
    make -j $(getconf _NPROCESSORS_ONLN) install

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
    echo "ldconfig" >> env.sh
