function [lId,l]= getLabove(lenArr,n) %% find index of subsequence length above 'n';lenArr are the retrieved array which contain all the subsequence length among 20 species

for i=1:length(lenArr)
templenArr=lenArr{i};
lId{i}=find(templenArr>n); %%lId: index
l{i}=templenArr(lId{i});   %% l: subsequence length
end

end


