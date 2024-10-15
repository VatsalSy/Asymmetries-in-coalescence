#!/bin/bash

Rr=("1e0" "1e0" "1e0" "1e0" "1e0" "1.5" "1.5" "1.5" "1.5" "1.5" "2.0" "2.0" "2.0" "2.0" "2.0" "4e0" "4e0" "4e0" "4e0" "4e0" "8e0" "8e0" "8e0" "8e0" "8e0")

start="1000"
end="1024"

for i in `seq $start $end`;
do

  echo $i
  scp -r *py $i/
  scp -r *get* $i/
  cd $i/
  
  python FacetsAndCOM.py $i ${Rr[$i-$start]}

  ffmpeg -framerate 300 -pattern_type glob -i 'VideoFacetsOnly/*.png' -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -c:v libx264 -r 30 -pix_fmt yuv420p $i-Facets.mp4 -y

  cd ../

done
