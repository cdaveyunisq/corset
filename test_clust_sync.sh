#!/bin/bash
mkdir -p examples/Cx_sitiens_sync
cp -f *.corset-recovery examples/Cx_sitiens_sync
cp -f build/corset_sync examples/Cx_sitiens_sync
pushd examples/Cx_sitiens_sync
echo "Running corset_sync in examples/Cx_sitiens_sync $(pwd)"
./corset_sync -f true -R true -p Cx_sitiens_sync ../Cx_sitiens_mosquitoes/W6.sorted.bam ../Cx_sitiens_mosquitoes/W9.sorted.bam
popd
