#!/bin/bash -e

if test $# -ne 1
then
  echo "usage: ./convertNoRelabel.sh <snap.txt>" 1>&2
  exit 1
fi

fname="${1}"

make

:<<END
if ! test -f "${fname}"
then
  for webgraph in ${dir}/*.properties
  do
    echo
    date
    printf "\nWEBGRAPH TO SNAP\n"
    java -cp $(echo ~/common/approxg/webgraph/*.jar|tr " " ":") it.unimi.dsi.webgraph.ArcListASCIIGraph -g BVGraph "${webgraph//.properties/}" "${fname}" || true
  done
fi
END

echo
date
printf "\nSNAP TO BIN\n"
bin/snapToBinNoRelabel "${fname}"
#rm "${fname}"

echo
date
printf "\nSORT\n"
bin/bsort -k 4 -r 8 "${fname}.rev.bin"

echo
date
printf "\nMAKE LISTS\n"
bin/makeLists "${fname}"
rm "${fname}.bin" "${fname}.rev.bin" "${fname}.raw.degree.bin"

echo
date
printf "\nCOMPACTIFY\n"
bin/compactify "${fname}"
rm "${fname}.raw.edge.bin" "${fname}.raw.vertex.bin"

