function [XfPlot,grp] = getHomoPlot2(tt,homoList,AAname) %%homology plot

for t=(1+tt*10):(10+tt*10)
    homoName=homoList{t};
    for an=1:18
    aaName=AAname{an};
    filename=[homoName,'homology',aaName,'.txt'];
    fileID=fopen(filename);
    Xhead=textscan(fileID,'%s %s %s %s %s\n',1,'Delimiter',',');
    XX=textscan(fileID,'%s %s %f %f %u\n','Delimiter',',');
    fclose(fileID);
    X{t,an}=XX{1,4};   %%% X: entropy location for 10 species within one group
    l(t,an)=length(X{t,an});
    end
end

    ll=sum(l);
    lll=cumsum(ll);
    XfPlot=cat(1,X{:});
    grp=strings([sum(ll),1]);
    grp(1:ll(1),:)=AAname{1};
    
    for i=2:18
    grp((1+lll(i-1)):lll(i),:)=AAname{i};
    end
        
end
