order={'11','12','13','14','15','21','22','23','31','32','33','34','35','36','37','41','42','43','44','45'};

AAlist={'Asp','Asn','Cys','Gln','Glu','His','Lys','Phe','Tyr','Ile','Ala','Gly','Pro','Thr','Val','Arg','Ser','Leu'};

SpeciesName={'Saccharomyces cerevisiae','Saccharomyces eubayanus','Saccharomyces kudriavzevii','Ashbya gossypii','Yarrowia lipolytica',...
    'Schizosaccharomyces pombe','Schizosaccharomyces octosporus','Schizosaccharomyces japonicus',...
    'Aspergillus clavatus','Aspergillus flavus','Aspergillus nidulans','Aspergillus niger','Aspergillus oryzae','Aspergillus terreus','Aspergillus fumigatus',...
    'Fusarium fujikori','Fusarium graminearum','Fusarium oxysporum','Fusarium verticilloides','Trichoderma reesei'};
% 
% for u=1:18
%     fileName2=([AAlist{u},'Accurate.txt']);
%     fileID=fopen(fileName2,'a');
%     fprintf(fileID,'%s,%s\n','SpeciesName','selectionPressure');
%     fclose(fileID);
%     
%     for o=1:20 %%%o number should be changed species counts
%         fileName1=([AAlist{u},'Accurate',order{o},'.mat']);
%         m=matfile(fileName1);
%         NHISTreal=m.NHISTreal;
%         SumHistF=m.SumHistF;
%         x=min(length(NHISTreal),length(SumHistF));
%         sumId=find((NHISTreal(1:x).*SumHistF(1:x))~=0);
%         selectP=sum(NHISTreal(sumId).*log(NHISTreal(sumId)./SumHistF(sumId))); %%%%Kullback Leibler divergence
%         
%         fileName2=([AAlist{u},'Accurate.txt']);
%         fileID=fopen(fileName2,'a');
%         fprintf(fileID,'%s,%d\n',SpeciesName{o},selectP);
%         fclose(fileID);
%     end
% end

X=zeros(20,18);
for u=1:18
%%  fileName3=([AAlist{u},'Accurate.txt']);
    fileName3=([AAlist{u},'.txt']);
    fileID=fopen(fileName3);
%     C0=textscan(fileID,'%s',2,'Delimiter',',');
    C=textscan(fileID,'%s %f','Delimiter',',');
    fclose(fileID);
    X(:,u)=C{1,2};
end

newO=[9,15,10,13,14,11,12,16,19,18,17,20,3,2,1,4,5,7,6,8];
XXp=X(newO,:);
save('KL20pre.mat','XXp');
figure
HeatMap(XXp,'RowLabels',SpeciesName(newO),'ColumnLabels',AAlist(:),'Colormap',[1 0.85 0.85;1 0.83 0.83;1 0.8 0.8;1 0.78 0.78;1 0.75 0.75;1 0.6 0.6;...
    1 0.55 0.55;1 0.4 0.4;1 0.35 0.35;1 0.31 0.31;1 0 0]);
figure
HeatMap(XXp,'RowLabels',SpeciesName(newO),'ColumnLabels',AAlist(:),'Colormap',[1 0.26 0.26;1 0.31 0.31;1 0.35 0.35;1 0.4 0.4;1 0.55 0.55;...
    1 0.6 0.6;1 0.65 0.65;1 0.7 0.7;1 0.75 0.75;1 0.86 0.86;1 0.9 0.9]);

% % % figure
% % % plot(edgesF,SumHistF,'g.');
% % % hold on
% % % plot(EdgesReal,NHISTreal,'r.');
% % % hold off