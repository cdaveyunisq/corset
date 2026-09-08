#!/bin/bash
mkdir -p examples/Cx_clean_u1_par
cp -f *.corset-recovery examples/Cx_clean_u1_par
cp -f build/corset_par examples/Cx_clean_u1_par
pushd examples/Cx_clean_u1_par
echo "Running corset_par in examples/Cx_clean_u1_par $(pwd)"
./corset_par -f true -R true -p Cx_clean_u1_par ../Cx_sitiens_mosquitoes/U1.clean.sorted.bam &>> corset_par.log
popd
