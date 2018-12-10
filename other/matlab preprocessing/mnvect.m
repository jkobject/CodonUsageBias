function FinPart = mnvect(Syno,SubLength) %% change the output of function 'partitions' to the format which is need for futher 'mnpdf' input
prepart=partitions(SubLength);
partResult = prepart(find(sum(prepart,2)<=Syno),:); %% find the total counts of each row which equals to required 'levels'
partResult = int32(partResult);
[Sy,Sx] = size(partResult);
FinPart = zeros(Sy,Syno);
for i=1:Sy
   for j=1:Sx
       if partResult(i,j)>0
           for k=1:partResult(i,j)
             FinPart(i,k+sum(partResult(i,1:j-1))) =j; %% change output foramt 
           end
       end
   end
   
end

        