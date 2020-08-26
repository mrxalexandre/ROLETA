#!/bin/bash


for i in log/*
do
	echo -n  "$i "
	grep solution $i | awk '{a+=$8} END {print a/NR}'
done

