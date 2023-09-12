#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
DATASET_OUTDIR=${SCRIPT_DIR}/TXT
echo "Download files into" ${DATASET_OUTDIR}

if test -f "${DATASET_OUTDIR}/wiki/snap.txt"; then
    echo "wiki exists"
else
    curl https://snap.stanford.edu/data/wiki-Vote.txt.gz --output wiki.txt.gz
    gzip -d wiki.txt.gz && mkdir -p ${DATASET_OUTDIR}/wiki && mv wiki.txt ${DATASET_OUTDIR}/wiki/snap.txt
fi

if test -f "${DATASET_OUTDIR}/youtube/snap.txt"; then
    echo "youtube exists"
else
    curl http://snap.stanford.edu/data/bigdata/communities/com-youtube.ungraph.txt.gz --output youtube.txt.gz
    gzip -d youtube.txt.gz && mkdir -p ${DATASET_OUTDIR}/youtube && mv youtube.txt ${DATASET_OUTDIR}/youtube/snap.txt
fi

if test -f "${DATASET_OUTDIR}/lj/snap.txt"; then
   echo "lj exists"
else
   curl http://snap.stanford.edu/data/bigdata/communities/com-lj.ungraph.txt.gz --output lj.txt.gz
   gzip -d lj.txt.gz && mkdir -p ${DATASET_OUTDIR}/lj && mv lj.txt ${DATASET_OUTDIR}/lj/snap.txt
fi

if test -f "${DATASET_OUTDIR}/orkut/snap.txt"; then
   echo "orkut exists"
else
   curl http://snap.stanford.edu/data/bigdata/communities/com-orkut.ungraph.txt.gz --output orkut.txt.gz
   gzip -d orkut.txt.gz && mkdir -p ${DATASET_OUTDIR}/orkut && mv orkut.txt ${DATASET_OUTDIR}/orkut/snap.txt
fi

if test -f "${DATASET_OUTDIR}/patents/snap.txt"; then
   echo "patents exitst"
else
   curl http://snap.stanford.edu/data/cit-Patents.txt.gz --output patents.txt.gz
   gzip -d patents.txt.gz && mkdir -p ${DATASET_OUTDIR}/patents && mv patents.txt ${DATASET_OUTDIR}/patents/snap.txt
fi

if test -f "${DATASET_OUTDIR}/friendster/snap.txt"; then
   echo "friendster exists"
else
   curl https://snap.stanford.edu/data/bigdata/communities/com-friendster.ungraph.txt.gz --output friendster.txt.gz
   gzip -d friendster.txt.gz && mkdir -p ${DATASET_OUTDIR}/friendster && mv friendster.txt ${DATASET_OUTDIR}/friendster/snap.txt
fi