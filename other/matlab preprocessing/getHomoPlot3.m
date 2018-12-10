function [XfPlot,grp] = getHomoPlot3(homoList,aaName) %%homology plot

for t=1:50
    homoName=homoList{t};
    filename=[homoName,'Rlhomology',aaName,'.txt'];
    fileID=fopen(filename);
    Xhead=textscan(fileID,'%s %s %s %s %s',1,'Delimiter',',');
    XX=textscan(fileID,'%s %s %f %f %u\n','Delimiter',',');
    fclose(fileID);
    X{t}=XX{1,4};    %%%%%% X: selection value
    l(t)=length(X{t});  %%%%% l: 50 homology counter  
end
    
    ll=cumsum(l);
    grp=strings([sum(l),1]);
    grp(1:l(1),:)=homoList{1}; %%% grp is for box plot
    
    for i=2:50
        grp((1+ll(i-1)):ll(i),:)=homoList{i};
    end
    
    XfPlot=cat(1,X{:});
end
