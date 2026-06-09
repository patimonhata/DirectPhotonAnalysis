#! /bin/bash
# export USER="$(id -u -n)"
# export LOGNAME=${USER}
# export HOME=/sphenix/u/${LOGNAME}

#source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.464
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
# source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.504


export MYINSTALL=/sphenix/u/ryotaro/DirectPhotonAnalysis/EmcalEtViewer/install

export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$MYINSTALL/include:$ROOT_INCLUDE_PATH

source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL


process=$1
nEvents=$2
save_tree=$3

echo process: ${process}
echo nEvents: ${nEvents}
echo save_tree: ${save_tree}

root.exe -q -b Fun4All_SingleParticlePi0.C\(${process},${nEvents},${save_tree}\)

echo all done process the Fun4All_DisplayLegoPlot.C, process_id: ${process}
