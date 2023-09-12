set -e

dir="${1}/"
fname="${dir}snap.txt"

make

./convert.sh ${dir}
./bin/labelize ${fname}
./bin/compactifylabel ${fname}
