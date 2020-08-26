#!/bin/bash

# dir com as instâncias
P=instances


# dir log
log="log";
mkdir $log;

for i in $P/*
do
	j=$(basename $i);

	echo "Running at $i";
	# Se o arquivo de log já existe, ele não executa
	if [ -e ${log}/${j}.log ];
	then
		echo "File exist..";
	else
		# Executa cada experimento 10x
		for seed in `seq 1 30`
		do
			# Aqui vc colocar a execução, jogando o a saida para o log
			# Adaptar aos parâmetros. No caso tempo -T com 30s e
			# a semente -s seed
			echo "Execução $seed " >> ${log}/${j}.log;
			./main -i $i -p 262 -e 0.1380 -m 0.3593 -r 0.6668 -K 5 -b 66 -x 4 -T 30 -s $seed >> ${log}/${j}.log;
			echo " " >> ${log}/${j}.log;

		done
	fi

done

