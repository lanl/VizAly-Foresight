# BLOSC
git clone https://github.com/Blosc/c-blosc.git
cd c-blosc/
git checkout v1.10.2
mkdir install
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..


# BigCrunch
git clone https://github.com/lanl/VizAly-BigCrunch.git
cd VizAly-BigCrunch/
mkdir install
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..


# SZ
git clone https://github.com/disheng222/SZ.git
cd SZ
mkdir install
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..