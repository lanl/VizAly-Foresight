echo "#!/bin/bash"

numNodes=$1
ranksPerNode=$2
projectPath=$3
inputPath=$4

echo "#!/bin/bash" > temp.sh
echo "#SBATCH -N $numNodes"         >> temp.sh
echo "#SBATCH --ntasks-per-node $ranksPerNode"  >> temp.sh
echo "#SBATCH -p scaling"           >> temp.sh
echo "#SBATCH -qos=scaling"         >> temp.sh

echo "source $projectPath/scripts/VizAly-CBench.bash.darwin"    >> temp.sh
echo "cd $projectPath/build"                                    >> temp.sh
echo "mpirun $projectPath/build/CBench $inputPath"              >> temp.sh

sbatch temp.sh
rm temp.sh

# Run as:
# darwinScalingJobCreator.sh 8 8 /projects/exasky/VizAly-CBench ../inputs/darwinStep499.json
# Param 1: num nodes
# Param 2: num ranks per node
# Param 3: project path
# Param 4: path to input
