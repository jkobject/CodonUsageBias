%% collect sample entropyMean 
parfor l=201:300

Pfin=sampEntropyMean6(6,l,10000000);

p1Sample6f201t300=log(Pfin);

end

save 'p1Sample6f201t300.mat' p1Sample6f201t300




parfor l=301:400

Pfin=sampEntropyMean6(6,l,10000000);

p1Sample6f301t400=log(Pfin);

end
save 'p1Sample6f301t400.mat' p1Sample6f301t400




parfor l=401:499

Pfin=sampEntropyMean6(6,l,10000000);

p1Sample6f401t499=log(Pfin);

end
save 'p1Sample6f401t499.mat' p1Sample6f401t499
