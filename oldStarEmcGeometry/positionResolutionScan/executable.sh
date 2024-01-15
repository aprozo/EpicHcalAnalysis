#!/bin/sh

export GUN_PHI=${3}
export GUN_THETA=${4}

export GUN_PARTICLE=neutron
export GUN_MOMENTUM=5
export NUMBER_OF_EVENTS=50000

export FILENAME=${GUN_PARTICLE}_${NUMBER_OF_EVENTS}events_p${GUN_MOMENTUM}gev_phi${GUN_PHI}_theta${GUN_THETA}
export DDSIM_FILE=sim_${FILENAME}.edm4hep.root
export EICRECON_FILE=eicrecon_${FILENAME}.edm4eic.root

echo "source /gpfs/mnt/gpfs02/eic/alpro/analysis/positionResolutionScan/runSimBatch.sh ${1} ${2}" | /eic/u/alpro/eic/eic-shell

#EOF
