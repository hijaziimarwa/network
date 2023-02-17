#!/bin/bash

seed=$1
num_fibers=$2

suffix="_seed"$seed"_numfibers"$num_fibers

paramfile=parameters$suffix
dir=current$suffix

rm $paramfile

cat parameters | while read line
do
    cont=$(echo $line | head -n1 | awk '{print $1;}')
    if [ "$cont" == "seed" ]
    then
        echo $cont $seed >> $paramfile
    elif [ "$cont" == "num_fibers" ]
    then
        echo $cont $num_fibers >> $paramfile
    else
	    echo $line >> $paramfile
    fi
done

mkdir ../$dir
cp main ../$dir/main
cp $paramfile ../$dir/$paramfile


cd ../$dir
./main $paramfile

cd ../
rm -rf dir
cd code_static

