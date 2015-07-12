#! /bin/bash

# usage: runcmake [-r|-d] [-n]
# default: release
# -r: release
# -d: debug
# -n: for note book

cd ./Script/

# release version
if   [ $# = 0 ]
then 
	cd ../Build/Release/
	cmake ../../
elif   [ $1 = '-r' ]
then
	cd ../Build/Release/
	if [ $# = 1 ]
	then
		cmake ../../ -DCMAKE_BUILD_TYPE=RELEASE
	elif [ $2='-n' ]
	then
		cmake ../../ -DCMAKE_BUILD_TYPE=MSSE
	fi
elif  [ $1 = '-n' ]
then
	cd ../Build/Release/
	cmake ../../ -DCMAKE_BUILD_TYPE=MSSE
	
# debug version
elif   [  $1 = '-d' ]
then
	cd ../Build/Debug/
	if [ $# = 1 ]
	then
		cmake ../../ -DCMAKE_BUILD_TYPE=DEBUG
	elif [ $2='-n' ]
	then
		cmake ../../ -DCMAKE_BUILD_TYPE=DMSSE
	fi
# report error
else
    echo "usage: runcmake [-r|-d] [-n]"
fi 