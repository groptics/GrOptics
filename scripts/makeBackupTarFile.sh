DATEM="28Sept3PM"

tar -cvf grOptics$DATEM.tar ./src ./include ./README ./scripts

gzip grOptics$DATEM.tar

scp grOptics$DATEM.tar.gz duke@ampere.math.grin.edu:~/GrOpticsBAK/.

rm grOptics$DATEM.tar.gz