# source environment
source env/bash.darwin.gcc-openmp

# move into the 3rdparty folder
root=$(eval pwd)
#cd source/3rdparty

# compile fpzip
#echo " "
#echo ===============================================================================================
#echo starting fpzip compilation
#echo ===============================================================================================
#echo " "
#cd fpzip-1.1.0/src
#make
#cd ../..

# compile sz
#echo " "
#echo ===============================================================================================
#echo starting sz compilation
#echo ===============================================================================================
#echo " "
#cd SZ-1.4.12.2
#make
#make install
#cd ..

# compile xz/lzma
#echo " "
#echo ===============================================================================================
#echo starting xz/lzma compilation
#echo ===============================================================================================
#echo " "
#cd xz-5.2.3
#make
#make install
#cd ..

# compile zfp
#echo " "
#echo ===============================================================================================
#echo starting zfp compilation
#echo ===============================================================================================
#echo " "
#cd zfp-0.5.0
#make
#cd ..

# create binary folder for cmake
echo " "
echo ===============================================================================================
echo starting cmake project compilation
echo ===============================================================================================
echo " "
cd $root/binary
make
make install
cd $root
