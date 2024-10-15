#!/bin/bash

start="1025"
end="1049"
for i in `seq $start $end`;
do
echo $i
mkdir -p $i
scp -r coal* $i/
scp -r *.h $i/
done
