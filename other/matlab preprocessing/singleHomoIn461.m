function X = singleHomoIn461(fileName,speciesName)

X=zeros(1,461);
X(:,:)=NaN;

fileID=fopen(fileName);
Xhead=textscan(fileID,'%s %s %s %s %s\n',1,'Delimiter',',');
XX=textscan(fileID,'%s %s %f %f %u\n','Delimiter',',');
fclose(fileID);

homoLen=length(XX{1,1});

HomoSpN=XX{1,1};
SeleVal=XX{1,4};

for l=1:homoLen
indWhole(l,1) = find(strcmp(speciesName, HomoSpN{l,1}));
end

X(indWhole)=SeleVal;

end