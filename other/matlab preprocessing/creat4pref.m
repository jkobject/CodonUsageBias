for l=101:300
if l<106
    tempM1=pf4f101t105{l};
    filename=['partition4o',num2str(l),'p.csv'];
    dlmwrite(filename,tempM1);
elseif l>105 && l<201
    tempM2=pf4f106t200{l};
    filename=['partition4o',num2str(l),'p.csv'];
    dlmwrite(filename,tempM2);
elseif l>200 && l<206
    tempM3=pf4f201t205{l};
    filename=['partition4o',num2str(l),'p.csv'];
    dlmwrite(filename,tempM3);
elseif l>205 && l<251
    tempM4=pf4f206t250{l};
    filename=['partition4o',num2str(l),'p.csv'];
    dlmwrite(filename,tempM4);
else
    tempM5=pf4f251t300{l};
    filename=['partition4o',num2str(l),'p.csv'];
    dlmwrite(filename,tempM5);
end
end
    
    
        