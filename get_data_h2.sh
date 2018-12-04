#!/bin/bash 

./get_bond_lengths_h2.sh > h2.bl 

grep "Total CI energy" h2_fci.out > h2.fci
grep "ERPA total energy" h2.out > h2.erp

paste h2.bl h2.fci h2.erp > h2.data 
awk 'BEGIN { OFS = "\t" } NR >= 1 { $12 = $11 - $6 } 1' h2.data> h2.data2
awk '{print $1,$12}' h2.data2>h2.data3
echo "Done collecting and sorting data ...."

