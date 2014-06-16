pwd
cd CMSSWDIR/src
pwd
eval `scramv1 runtime -sh`
cd -
pwd
cd BASEDIR
pwd
source scripts/setup.sh
cd JOBDIR
pwd
unbuffer EXENAME   studyNoise_LABEL.cfg > JOBDIR/out.txt
