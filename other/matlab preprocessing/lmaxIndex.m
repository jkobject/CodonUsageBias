function [lmaxId,lmax] = lmaxIndex(lenArr)

for i=1:length(lenArr)
lenArrTemp=lenArr{i};
MaxId=find(lenArrTemp==(max(lenArrTemp)));
lmaxId(i)=MaxId(1);
lmax(i)=lenArrTemp(MaxId(1));
end
end