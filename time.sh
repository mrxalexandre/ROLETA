#!/bin/bash


for i in log/*
do
	j="${i/log/} "
	j="${j/\//} "
	j="${j/.dat.log/} "
	echo -n "$j "
	grep Total $i | awk '{a+=$4} END {printf a/NR}'
	echo  -n " "
	grep ttb $i | awk '{a+=$6} END {print a/NR}'
done

