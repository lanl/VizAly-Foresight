#!/bin/bash

rootDir=`pwd`

NODE_URL=https://nodejs.org/dist/v6.11.0/node-v6.11.0.tar.gz
NODE_SRC_DIR=$rootDir/node-v6.11.0
NODE_DIR=$rootDir/node-v6.11.0-install

# download and install node.js
cd $rootDir
if [ ! -d "$NODE_DIR" ] ; then
  curl -L $NODE_URL | tar zxf -
  cd $NODE_SRC_DIR
  ./configure --prefix=$NODE_DIR
  make
  make install
fi


# download z-checker-web
if [ ! -d z-checker-web ] ; then
  git clone https://github.com/CODARcode/z-checker-web
else
  cd z-checker-web
  git pull
  cd ..
fi

export PATH=$NODE_DIR/bin:$PATH
cd z-checker-web
npm install
npm install -g http-server

cd $rootDir
