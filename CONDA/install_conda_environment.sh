#!/bin/bash
## conda environment name

source /data/shamsaddinisha/conda/etc/profile.d/conda.sh

ENV_NAME_PY2=ATAC_Seq_py2
ENV_NAME_PY3=ATAC_Seq_py3

PY2_Version=python=2.7.15

PY3_Version=python=3.7.0


if [[ $(conda env list | grep $ENV_NAME_PY2 | wc -l) == "0" ]]; then
	echo "Installing $ENV_NAME_PY2 ...."
	conda create -n $ENV_NAME_PY2 $PY2_Version -y -c defaults
	echo "$ENV_NAME_PY2 Installed!!!"
else
	echo "$ENV_NAME_PY2 is Available"
fi


if [[ $(conda env list | grep $ENV_NAME_PY3 | wc -l) == "0" ]]; then
	echo "Installing $ENV_NAME_PY3 ...."
	conda create -n $ENV_NAME_PY3 $PY3_Version -y -c defaults
	echo "$ENV_NAME_PY3 Installed!!!"
else
	echo "$ENV_NAME_PY3 is Available"
fi


echo "Installing environment has been successfully done."