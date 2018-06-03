#!/bin/sh
for i in *.c
do
 echo $i 
 diff $i $1$i
done
