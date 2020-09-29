#! /bin/bash

PYTHONPATH="/localhome/py13je/work/athena-cwb-dust/vis/python"

oldfname=""

while true
do
  fname=`ls -Art *athdf | tail -n 1`
  if [ "$fname" == "$oldfname" ]; then
    echo "Writing new image from dataset $fname"
    python $PYTHONPATH/plot_slice.py $fname dens test.png --logc
    echo "Done!"
  fi
  oldfname=`echo $fname`
done
