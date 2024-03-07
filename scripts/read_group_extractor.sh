#!/bin/bash

#set -ueo pipefail
#set -x

header=$(zcat $1 | head -n 1)
echo $header
id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
echo "id $id"
sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+$")
echo "sm  $sm"
rg="./@RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA"

echo $rg