%% in one species for all the 18 amino acids after replacement

for i=1:18

xa{i}= MidFq1{1,i}{:,1};
ya{i}= MidFq1{1,i}{:,2};

% xaT{i}= MidFqT6{1,i}{:,1};
% yaT{i}= MidFqT6{1,i}{:,2};

end
for i=1:18
plot(xb{i},xa{i},'o')
end