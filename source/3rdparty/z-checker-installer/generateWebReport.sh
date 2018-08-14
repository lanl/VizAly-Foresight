#!/bin/bash

if [[ $# < 1 ]]
then
	echo "Usage: $0 <casename>"
	exit
fi 

testcase=$1

rootDir=`pwd`

export PATH=$rootDir/gnuplot-5.0.6-install/bin:$rootDir/node-v6.11.0-install/bin:$PATH
# export PATH=$rootDir/node-v6.11.0-install/bin:$PATH

cd $rootDir/z-checker-web
# npm update
node generate_report.js $rootDir/Z-checker/$testcase

echo "Please use your web browser to open this webpage: http://localhost:8087"
echo "Press Ctrl+C to exit."

cd $rootDir/Z-checker/$testcase/report/web
http-server -p8087

# echo "file://"$rootDir/Z-checker/$testcase/report/web/index.html
