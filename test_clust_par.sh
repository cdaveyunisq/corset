#!/bin/bash
mkdir -p examples/Cx_sitiens_par
cp -f *.corset-recovery examples/Cx_sitiens_par
cp -f build/corset_par examples/Cx_sitiens_par
pushd examples/Cx_sitiens_par
echo "Running corset_par in examples/Cx_sitiens_par $(pwd)"
./corset_par -f true -R true -p Cx_sitiens_par ../Cx_sitiens_mosquitoes/W6.sorted.bam ../Cx_sitiens_mosquitoes/W9.sorted.bam
popd
