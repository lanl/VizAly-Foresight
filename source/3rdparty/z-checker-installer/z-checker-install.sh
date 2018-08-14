#!/bin/bash

rootDir=`pwd`

#----------check X11------------------------
X_PATH=`which X`
if [ ! -x ${X_PATH} ];then
	echo "Error: missing X11!"
	echo "Please install X11 first (requiring root previlege)."
	exit
fi

#----------check git and perl---------------
PERL_PATH=`which perl`
if [ ! -x ${PERL_PATH} ]; then
	echo "Error: missing perl command; Please install perl."
	exit
fi

#---------- download gnuplot ----------------
GNUPLOT_URL="https://downloads.sourceforge.net/project/gnuplot/gnuplot/5.0.6/gnuplot-5.0.6.tar.gz"
GNUPLOT_SRC_DIR=$rootDir/gnuplot-5.0.6
GNUPLOT_DIR=$rootDir/gnuplot-5.0.6-install

GNUPLOT_EXE_PATH=`which gnuplot`
if [ ! -x "$GNUPLOT_EXE_PATH" ]; then
	if [ ! -d "$GNUPLOT_DIR" ]; then
		# download gnuplot source
		curl -L $GNUPLOT_URL | tar zxf -
		if [ ! -d "$GNUPLOT_SRC_DIR" ] ; then
			echo "FATAL: cannot download and extract gnuplot source."
			exit
		fi

		# compile gnuplot
		cd $GNUPLOT_SRC_DIR
		./configure --prefix=$GNUPLOT_DIR
		make && make install
		cd $rootDir
		echo "export GNUPLOT_HOME=$GNUPLOT_DIR" > $rootDir/env_config.sh
		echo "export PATH=\$PATH:\$GNUPLOT_HOME/bin" >> $rootDir/env_config.sh
	fi

fi

#---------- download Z-checker --------------
cd $rootDir
git clone https://github.com/CODARcode/Z-checker.git
cd Z-checker
./configure --prefix=$rootDir/Z-checker/zc-install
make
make install
cp ../zc-patches/generateReport.sh ./examples/

#---------- download ZFP and set the configuration -----------
cd $rootDir

git clone https://github.com/LLNL/zfp.git
cd zfp
make

cd -
cp zfp-patches/zfp-zc.c zfp/utils
cp zfp-patches/*.sh zfp/utils

cd zfp/utils/
patch -p0 < ../../zfp-patches/Makefile-zc.patch
make

#---------- download SZ and set the configuration -----------
cd $rootDir
git clone https://github.com/disheng222/SZ

cd SZ/sz/src
patch -p1 < ../../../sz-patches/sz-src-hacc.patch

cd ../..
./configure --prefix=$rootDir/SZ/sz-install
make
make install

cd example
patch -p0 < ../../sz-patches/Makefile-zc.bk.patch
make -f Makefile.bk
cp ../../Z-checker/examples/zc.config .
patch -p0 < ../../zc-patches/zc-probe.config.patch

cp ../../sz-patches/sz-zc-ratedistortion.sh .
cp ../../sz-patches/testfloat_CompDecomp.sh .
cp ../../sz-patches/testdouble_CompDecomp.sh .

#----------- download latexmk --------------------------------
cd $rootDir
latexmk_url=http://ctan.math.utah.edu/ctan/tex-archive/support/latexmk.zip
latexmk_dir=latexmk
latexmk_exe_path=`which latexmk`
if [ ! -x "$latexmk_exe_path" ]; then
	if [ ! -d "$latexmk_dir" ]; then
		curl -O $latexmk_url
		unzip latexmk.zip
		cd $latexmk_dir
		ln -s "$rootDir/$latexmk_dir/latexmk.pl" latexmk
		echo "export LATEXMK_HOME=$rootDir/$latexmk_dir" >> $rootDir/env_config.sh
		echo "export PATH=\$PATH:\$LATEXMK_HOME" >> $rootDir/env_config.sh
		cd $rootDir
		rm -rf latexmk.zip
	fi
fi

#----------- download ghost view (gsview) if necessary-----------
cd $rootDir
ghost_url="https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs921/ghostpdl-9.21.tar.gz"
ghost_src_dir=$rootDir/ghostpdl-9.21
ghost_dir=$rootDir/ghostpdl-9.21-install
PS2PDF_EXE_PATH=`which ps2pdf`
if [ ! -x "$PS2PDF_EXE_PATH" ]; then
        if [ ! -d "$ghost_dir" ]; then
                # download ghost source
                curl -L $ghost_url | tar zxf -
                if [ ! -d "$ghost_src_dir" ] ; then
                        echo "FATAL: cannot download and extract ghost source."
                        exit
                fi

                # compile ghost
                cd $ghost_src_dir
                ./configure --prefix=$ghost_dir
                make && make install
                cd $rootDir
		echo "export GHOST_HOME=$ghost_dir" >> $rootDir/env_config.sh
		echo "export PATH=\$PATH:\$GHOST_HOME/bin" >> $rootDir/env_config.sh
        fi

fi


if [ -f $rootDir/env_config.sh ]; then
	mv $rootDir/env_config.sh $rootDir/Z-checker/examples/env_config.sh
fi
