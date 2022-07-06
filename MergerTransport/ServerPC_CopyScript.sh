#!/bin/bash

if [[ $# -eq 0 ]]
then
	echo "usage: ./CopyScript FolderToCopy"
	exit 1
fi

if [[ $# -eq 1 ]]
then
    para=$1
fi

pathMath="/users/waleed.khalid/Mma/MergerTransport/$para"
scp -r waleed.khalid@cbe.vbc.ac.at:$pathMath .


