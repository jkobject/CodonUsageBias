clc 
clear
%% distribution plot for Gly (current sequence from input & with replacement)

ctG={'GGG','GGA','GGT','GGC'}; %%Gly
%cfG=[0.12,0.23,0.45,0.20];
cfG=[1/4,1/4,1/4,1/4];

pasteCodon=getCodonSequence; %% call function getCodonSequence to creat input 


for i=1:length(pasteCodon)
    
    [Ng,G(i,:)]=GlyH(pasteCodon{1,i});
    
    r = rand(1,Ng);

for j=1:Ng
    if r(j)<cfG(1)
        SgR(j)=ctG(1);
        
    elseif r(j)>=cfG(1) && r(j)<(cfG(1)+cfG(2))
        SgR(j)=ctG(2);
        
    elseif r(j)>=(cfG(1)+cfG(2)) && r(j)<(cfG(1)+cfG(2)+cfG(3))
         SgR(j)=ctG(3);
    else SgR(j)=ctG(4);
    end
end

[NRg,GR(i,:)]=GlyH((SgR)');

end
   


