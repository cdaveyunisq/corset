#!/bin/bash
mkdir -p examples/Cx_sitiens_par1
cp -f *.corset-recovery examples/Cx_sitiens_par1
cp -f build/corset_par examples/Cx_sitiens_par1
pushd examples/Cx_sitiens_par1
echo "Running corset_par in examples/Cx_sitiens_par1 $(pwd)"
./corset_par -f true -R true -p Cx_sitiens_par1 ../Cx_sitiens_mosquitoes/W6.sorted.bam ../Cx_sitiens_mosquitoes/W9.sorted.bam &>> corset_par.log
popd
