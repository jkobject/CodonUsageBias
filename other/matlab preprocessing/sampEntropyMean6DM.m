function [] = sampEntropyMean6DM(Syno,Length,iter) %% Syno: number of synonymous codons; Length: sequence length; iter:interation counts
P=zeros(1,Syno);
P(1:Syno)=1/Syno;
initlX=zeros(1,Syno);  %% initlX : matrix of partitions
initlX(1)=Length;
% initlX=[7,7,6];

filename=['Tst',num2str(Length),'p.csv'];
fileID = fopen(filename,'w');
for i=1:iter
OutX1=initlX;
p1=mnpdf(initlX,P);    
    rFm=rand;  %% rFm: random number to decide which codon to minus 1;
    s=randi(6); %% sFa: random number to decide which codon to add 1;
    
    r1=initlX(1)/Length;
    r2=(initlX(1)+initlX(2))/Length;
    r3=(initlX(1)+initlX(2)+initlX(3))/Length;
    r4=(initlX(1)+initlX(2)+initlX(3)+initlX(4))/Length;
    r5=(initlX(1)+initlX(2)+initlX(3)+initlX(4)+initlX(5))/Length;
    
    if rFm<=r1
        t=1;
    elseif (r1<rFm)&&(rFm<=r2)
            t=2;
    elseif (r2<rFm)&&(rFm<=r3)
            t=3;
    elseif (r3<rFm)&&(rFm<=r4)
            t=4;
    elseif (r4<rFm)&&(rFm<=r5)
            t=5;        
    else    t=6;
   % else    t=3;
            
    end
    
    
    initlX(t)=initlX(t)-1;
    
    
    
    %s=sFa;   %% t: codon index to minus; s: codon index to add
            
    initlX(s)=initlX(s)+1;
    
    
    p2=mnpdf(initlX,P);
%%    p(i)=mnpdf(initlX,P);          %% p(i) store sample mnpdf values, ith is the current one
   Acp=min((p2.*((initlX(t)+1)/Length))/(p1.*(initlX(s)/Length)),1); 
   
   
   %  Acp(i)=min((p(i).*((initlX(t)+1)/Length).*(1/Syno))/(p(i-1).*(initlX(s)/Length).*(1/Syno)),1); %%Acp(i): acceptance probability, if > 1 accept
    
    if Acp > rand %% whether or not accept rand(1)
        Pfin=p2;
        OutX2=initlX;
 %      disp(Acp); disp("\n");
  %     disp(initlX);
    else
        Pfin=p1;
        OutX2=OutX1;
        initlX=OutX1;
    end

%fprintf(fileID,'%d,',Pfin); 
 fprintf(fileID,'%d',Pfin); 

 fprintf(fileID,',%d',initlX); 
   fprintf(fileID,'\n'); 
end
fclose(fileID); 
% Xwrite=cell2mat(OutX);   %% use the following code to store P data in Pwrite.txt
% fileID = fopen('Xwrite.txt','w');
% % fprintf(fileID,'%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f\n',Xwrite);
% fprintf(fileID,'%3i,%3i,%3i,%3i\n',Xwrite);
% fclose(fileID);

% figure
% histogram(cell2mat(Pfin));

end
