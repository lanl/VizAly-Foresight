# source environment
source env/bash.darwin.gcc-openmp

# move into the 3rdparty folder
root=$(eval pwd)
cd source/3rdparty

# clean fpzip
echo " "
echo ===============================================================================================
echo starting fpzip cleanup
echo ===============================================================================================
echo " "
cd fpzip-1.1.0/src
make clean
cd ../..

# clean sz
echo " "
echo ===============================================================================================
echo starting sz cleanup
echo ===============================================================================================
echo " "
cd sz-1.4.9-hacc
make clean
make distclean
cd ..

# clean xz/lzma
echo " "
echo ===============================================================================================
echo starting xz/lzma cleanup
echo ===============================================================================================
echo " "
cd xz-5.2.3
make clean
make distclean
cd ..

# clean zfp
echo " "
echo ===============================================================================================
echo starting zfp cleanup
echo ===============================================================================================
echo " "
cd zfp-0.5.0
make clean
cd ..

# clean cmake project
echo " "
echo ===============================================================================================
echo starting cmake project cleanup
echo ===============================================================================================
echo " "
cd $root
cd binary
make clean
cd $root

if [ -d binary ]
then
	rm -rf binary
fi

if [ -d build ]
then
	rm -rf build
fi
