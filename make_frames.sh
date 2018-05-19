#!/bin/bash


fn=$1


for count in {0..9}
do
   convert ${fn}000${count}.ppm frame_000${count}.jpg
done

for count in {10..99}
do
   convert ${fn}00${count}.ppm frame_00${count}.jpg
done

