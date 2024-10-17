#!/bin/bash

# check code reproduces results on two months of jb1 data
path='../src/aodpy/'
python ${path}envcal.py -t

python ${path}langley.py -t

python ${path}gencal.py -t

python ${path}aod2p.py -t

if diff -qr SampleData CorrectOutput > /dev/null; then
	echo "Test Passed"
else
        diff -qr SampleData CorrectOutput > testresults.txt
	echo "Test Failed, see testresults.txt"
fi
