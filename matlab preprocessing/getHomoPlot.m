function [XfPlot,grp] = getHomoPlot(tt,homoList,aaName) %%homology plot

for t=(1+tt*10):(10+tt*10)
    homoName=homoList{t};
    filename=[homoName,'homology',aaName,'.txt'];
    fileID=fopen(filename);
    Xhead=textscan(fileID,'%s %s %s %s %s',1,'Delimiter',',');
    XX=textscan(fileID,'%s %s %f %f %u\n','Delimiter',',');
    fclose(fileID);
    X{t-tt*10}=XX{1,4};   %%% X: entropy location for 10 species within one group
    l(t-tt*10)=length(X{t-tt*10});  %%%% l:counts for each of 10 species within one group
end
    
    ll=cumsum(l);
    grp=strings([sum(l),1]);
    grp(1:l(1),:)=homoList{1+tt*10}; %%% grp is for box plot
    
    for i=2:10
        grp((1+ll(i-1)):ll(i),:)=homoList{i+tt*10};
    end
    
    XfPlot=cat(1,X{:});
end
