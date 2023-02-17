#!/bin/bash

paramfile=$1
dir=$2

mkdir ../$dir
cp main ../$dir/main
cp $paramfile ../$dir/$paramfile


cd ../$dir
./main $paramfile

cd ../
rm -rf dir
cd code_static

