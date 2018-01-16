for l=1:400
if l<201
    tempM1=pf2f1t200{l};
    filename=['partition2o',num2str(l),'p.csv'];
    dlmwrite(filename,tempM1);
else
    tempM2=pf2f201t400{l};
    filename=['partition2o',num2str(l),'p.csv'];
    dlmwrite(filename,tempM2);
end
end
    
    
        