parfor l=101:200       %% myrtle calculate p : 6 synonymous codons of length 101 to 200

Pfin=sampSequenceD(6,l,10000);

p1Sample6{l}=log(Pfin);

end

save 'p1Sample6f101t200.mat' p1Sample6