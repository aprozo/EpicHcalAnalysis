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

if [ -f "$DDSIM_FILE" ]; then
    echo "$DDSIM_FILE exists."
else
    source /gpfs/mnt/gpfs02/eic/alpro/epic/install/setup.sh
    ddsim --compactFile $DETECTOR_PATH/epic.xml --numberOfEvents ${NUMBER_OF_EVENTS} --random.seed $(date +%N) --enableGun --gun.particle neutron --gun.thetaMin ${GUN_THETA_MIN}*degree --gun.thetaMax ${GUN_THETA_MAX}*degree --gun.phiMin ${GUN_PHI_MIN}*degree --gun.phiMax ${GUN_PHI_MAX}*degree --gun.distribution uniform --gun.momentumMin ${GUN_MOMENTUM_MIN}*GeV --gun.momentumMax ${GUN_MOMENTUM_MAX}*GeV --outputFile ${DDSIM_FILE}
fi

if [ -f "$EICRECON_FILE" ]; then
    echo "$EICRECON_FILE exists."
else
    source /gpfs/mnt/gpfs02/eic/alpro/EICrecon/bin/eicrecon-this.sh
    eicrecon $DDSIM_FILE -Ppodio:output_file=${EICRECON_FILE} -Pjana:nevents=${NUMBER_OF_EVENTS} -Ppodio:output_include_collections="
    MCParticles,HcalEndcapNRawHits,HcalEndcapNRecHits,HcalEndcapNMergedHits,HcalEndcapNClusters,HcalEndcapNTruthClusters,EcalEndcapNRawHits,
    EcalEndcapNRecHits,EcalEndcapNClusters,EcalEndcapNClusterAssociations,EcalEndcapNTruthClusters,EcalEndcapNTruthClusterAssociations, "
    rm -r calibrations
    rm -r fieldmaps
fi

root -l '/gpfs/mnt/gpfs02/eic/alpro/analysis/positionResolutionScan/readHCalRecoReader.C("'${EICRECON_FILE}'" , "outputFill.root")'
