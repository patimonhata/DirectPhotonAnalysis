#! /bin/bash
# export USER="$(id -u -n)"
# export LOGNAME=${USER}
# export HOME=/sphenix/u/${LOGNAME}

#source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.464
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
# source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.532


export MYINSTALL=/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install

export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$MYINSTALL/include:$ROOT_INCLUDE_PATH

source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL


process=$1
nEvents=$2

echo process: ${process}
echo nEvents: ${nEvents}

root.exe -q -b Fun4All_Pi0Reconstruction.C\(${process},${nEvents}\)

echo all done process the Fun4All_Pi0Reconstruction.C, process_id: ${process}
