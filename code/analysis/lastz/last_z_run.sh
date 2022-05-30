#!/usr/bin/bash

sbatch lastz.sh \
../data/ortholog/chicken_assembly/chicken_repeat/fix_sequences.fa.masked \
../../../data/genome/repeatmask/Struthio_camelus.20130116.OM.fa.masked \
../../../data/lastz/chicken_ostrich.lastz
