#!/bin/sh
#
#SBATCH --job-name="RUN_CREST"
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1GB
#SBATCH --account=education-as-msc-ce

# Collect parsed variables through python

path_to_complexes=$1
method=$2
charge=$3
geom=$4
multiplicity=$5
solvent=$6

cd "${path_to_complexes}"

FILES=`ls *.xyz`

for file in $FILES; do
filename=${file::-4}
mkdir -p ${filename}
cd $filename
# copy the org xyz to file var
cd .. 
cp $file $filename
cd $filename

if [[ "${geom}" == *"SP"* ]]; then

echo "charge_${file}=${charge}"
echo "Running a pre-opt with XTB"
xtb ${file} --${method} --chrg $charge --uhf 1 --alpb ${solvent} --opt 1 > ${file}_xtb.out
echo "Run crest from xtbopt.xyz"
# crest xtbopt.xyz --gfn2 --T 8 --chrg $charge --uhf 1 --alpb ${solvent} > CREST.out
else
echo "charge_${file}=${charge}"
echo "Running a pre-opt with XTB"
xtb ${file} --${method} --chrg $charge --uhf ${multiplicity} --alpb ${solvent} --opt 1 > ${filename}_xtb.out
echo "Run crest from pre-opt"
# crest xtbopt.xyz --${method} --T 8 --chrg $charge --uhf ${multiplicity} --alpb ${solvent} > CREST.out
fi
echo "CREST done"
cd ..
mv $filename ../CREST/${filename}
done