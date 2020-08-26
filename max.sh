#!/bin/bash


for i in log/*
do
	echo -n  "$i "
	grep solution $i | awk 'BEGIN { a=-1} {if(a<$9) a=$9} END {print a}'
done

