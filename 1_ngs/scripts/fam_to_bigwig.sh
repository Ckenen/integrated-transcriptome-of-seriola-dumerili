#!/bin/bash

# Usage: famToBigwig.sh --fr input.fam prefix
strand=$1
fam=$2
prefix=$3

if [ $strand == "--fr" ]; then
    forward="+"
    reverse="-"
elif [ $strand == "--rf" ]; then
    forward="-"
    reverse="+"
else
    echo "error! strand is --fr or --rf"
    exit 1
fi

# GENOME SIZES
gsize="${prefix}.famToBigwig.TEMP.sizes"
samtools view -H ${fam} | grep '@SQ' | awk -v OFS='\t' '{print $2,$3}' | sed 's/SN://g' | sed 's/LN://g' | sort -k1,1 > ${gsize}

# LIBRARY SIZE
bed="${prefix}.famToBigwig.TEMP.bed"
famtools.py convert -f bed $fam $bed
lsize=`wc -l $bed | awk '{print $1}'`
scale=`echo ${lsize} | awk '{print 1000000/$0}'`
# echo "Library size: "$lsize
# echo "Scale: "$scale


# RAW
bg="${prefix}.famToBigwig.TEMP.bedGraph"
bw="${prefix}.bw"
bedtools genomecov -bg -split -i ${bed} -g ${gsize} | sort -k1,1 -k2,2n > ${bg}
bedGraphToBigWig ${bg} ${gsize} ${bw}

bg="${prefix}.famToBigwig.TEMP.bedGraph"
bw="${prefix}.+.bw"
bedtools genomecov -bg -split -strand ${forward} -i ${bed} -g ${gsize} | sort -k1,1 -k2,2n > ${bg}
bedGraphToBigWig ${bg} ${gsize} ${bw}

bg="${prefix}.famToBigwig.TEMP.bedGraph"
bw="${prefix}.-.bw"
bedtools genomecov -bg -split -strand ${reverse} -i ${bed} -g ${gsize} | sort -k1,1 -k2,2n > ${bg}
bedGraphToBigWig ${bg} ${gsize} ${bw}


# NORMALIZED
bg="${prefix}.famToBigwig.TEMP.bedGraph"
bw="${prefix}.normalized.bw"
bedtools genomecov -bg -split -scale ${scale} -i ${bed} -g ${gsize} | sort -k1,1 -k2,2n > ${bg}
bedGraphToBigWig ${bg} ${gsize} ${bw}

bg="${prefix}.famToBigwig.TEMP.bedGraph"
bw="${prefix}.normalized.+.bw"
bedtools genomecov -bg -split -scale ${scale} -strand ${forward} -i ${bed} -g ${gsize} | sort -k1,1 -k2,2n > ${bg}
bedGraphToBigWig ${bg} ${gsize} ${bw}

bg="${prefix}.famToBigwig.TEMP.bedGraph"
bw="${prefix}.normalized.-.bw"
bedtools genomecov -bg -split -scale ${scale} -strand ${reverse} -i ${bed} -g ${gsize} | sort -k1,1 -k2,2n > ${bg}
bedGraphToBigWig ${bg} ${gsize} ${bw}


# CLEAN
rm $gsize $bed $bg

# echo "Finished!"

