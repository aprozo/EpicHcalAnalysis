#!/bin/sh

workdir=/gpfs/mnt/gpfs02/eic/alpro/output/positionResolutionScan/Phi${GUN_PHI}Theta${GUN_THETA}
mkdir -p $workdir
cd $workdir

export GUN_THETA_MIN=$(echo "$GUN_THETA - 0.0001" | bc)
export GUN_THETA_MAX=$(echo "$GUN_THETA + 0.0001" | bc)
export GUN_PHI_MIN=$(echo "$GUN_PHI - 0.0001" | bc)
export GUN_PHI_MAX=$(echo "$GUN_PHI + 0.0001" | bc)
export GUN_MOMENTUM_MIN=$(echo "$GUN_MOMENTUM - 0.00001" | bc)
export GUN_MOMENTUM_MAX=$(echo "$GUN_MOMENTUM + 0.00001" | bc)

source /gpfs/mnt/gpfs02/eic/alpro/epic/install/setup.sh
ddsim --compactFile $DETECTOR_PATH/epic.xml --numberOfEvents ${NUMBER_OF_EVENTS} --random.seed $(date +%N) --enableGun --gun.particle neutron --gun.thetaMin ${GUN_THETA_MIN}*degree --gun.thetaMax ${GUN_THETA_MAX}*degree --gun.phiMin ${GUN_PHI_MIN}*degree --gun.phiMax ${GUN_PHI_MAX}*degree --gun.distribution uniform --gun.momentumMin ${GUN_MOMENTUM_MIN}*GeV --gun.momentumMax ${GUN_MOMENTUM_MAX}*GeV --outputFile ${DDSIM_FILE}

source /gpfs/mnt/gpfs02/eic/alpro/EICrecon/bin/eicrecon-this.sh
eicrecon $DDSIM_FILE -Ppodio:output_file=${EICRECON_FILE} -Pjana:nevents=${NUMBER_OF_EVENTS} -Ppodio:output_include_collections="MCParticles,HcalEndcapNRawHits,HcalEndcapNRecHits,HcalEndcapNMergedHits,HcalEndcapNClusters,HcalEndcapNTruthClusters,EcalEndcapNRawHits,EcalEndcapNRecHits,EcalEndcapNClusters,EcalEndcapNClusterAssociations,EcalEndcapNTruthClusters,EcalEndcapNTruthClusterAssociations"

root -l '/gpfs/mnt/gpfs02/eic/alpro/analysis/positionResolutionScan/readHCalRecoReader.C("'${EICRECON_FILE}'" , "outputFill.root")'
root -l '/gpfs/mnt/gpfs02/eic/alpro/analysis/positionResolutionScan/fithistos_1D.C("outputFill.root")'

rm -r calibrations
rm -r fieldmaps