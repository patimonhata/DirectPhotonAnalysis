#! /bin/bash
# export USER="$(id -u -n)"
# export LOGNAME=${USER}
#source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.464
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
# source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.532


export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$MYINSTALL/include:$ROOT_INCLUDE_PATH

# source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL


process=$1

echo process: ${process}

root.exe -q -b MakeTruthPi0HistogramsFromEventDisplayTree.C\(${process}\)

echo all done process the MakeTruthPi0HistogramsFromEventDisplayTree.C, process_id: ${process}
