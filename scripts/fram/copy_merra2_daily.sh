#!/bin/bash

dir1=/home/jwalker/eady/datastore/merra2/daily
dir2=/home/jwalker/datastore/merra2/daily

for yr in {1980..2015}
do
    echo `date`
    echo $yr
    src=$dir1/$yr
    dest=$dir2
    echo Copying $src to $dest
    cp -r $src $dest
done
echo 'Done!'
