for l=1:400
if l<201
    tempM1=pf3f1t200{l};
    filename=['partition3o',num2str(l),'p.csv'];
    dlmwrite(filename,tempM1);
else
    tempM2=pf3f201t400{l};
    filename=['partition3o',num2str(l),'p.csv'];
    dlmwrite(filename,tempM2);
end
end
    
    
        