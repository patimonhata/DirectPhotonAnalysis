#! /bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana

export MYINSTALL=/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install
export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$MYINSTALL/include:$ROOT_INCLUDE_PATH

source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

process=$1
nEvents=$2

echo process: ${process}
echo nEvents: ${nEvents}

root.exe -q -b Fun4All_Pi0ReconstructionSplitComparison.C\(${process},${nEvents}\)
status=$?

if [ ${status} -ne 0 ]; then
  echo Fun4All_Pi0ReconstructionSplitComparison.C failed with status ${status}
  exit ${status}
fi

echo all done process the split comparison, process_id: ${process}
