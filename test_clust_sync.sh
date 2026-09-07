#!/bin/bash
mkdir -p examples/Cx_sitiens_sync
build/corset_sync -f true -R true -p examples/Cx_sitiens_sync examples/Cx_sitiens_mosquitoes/W6.sorted.bam
