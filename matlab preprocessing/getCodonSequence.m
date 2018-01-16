    function pasteCodonArray = getCodonSequence(filename)

    %%read the sequence data from file 'inputSequence'
%  
%     fileID = fopen(filename);
%     tt = textscan(fileID,'%s','delimiter',';');
%     fclose(fileID);
% 
%     %%seperate the sequence data name & sequence itself
%     xx = tt{1};
%     l = length(xx);
%     datRaw = xx(2:2:l,:);   
% 
%     for i=1:length(datRaw)
% 
%         datPre=datRaw{i};
% 
%         if mod(length(datPre),3)~=0
%             
%         k=strfind(datPre,'ATG');  %% k(1): the indice of first start codon
% 
%         if numel(k)~=0    %% saeu,saku,sj; have the case numel=0, why??  
% 
%         datPre=datPre(k(1):length(datPre));
% 
%         u=strfind(datPre,'TAA'); %% m: the indice of first stop codon
% 
%         v=strfind(datPre,'TGA');
% 
%         w=strfind(datPre,'TAG');
% 
%         M=[u,v,w];
% 
%         m=sort(M);
% 
%         n=find (~(mod((m+2),3)),1); %% n: indice in m, choose the first proper stop codon which can generate the sequence length as 3 multiples
%           
%         if n==0
%               
%               datPre2=datPre(k(2):length(datPre));
% 
%               u2=strfind(datPre2,'TAA'); %% m: the indice of first stop codon
% 
%               v2=strfind(datPre2,'TGA');
% 
%               w2=strfind(datPre2,'TAG');
% 
%               M2=[u2,v2,w2];
% 
%               m2=sort(M2);
%               
%               n2=find (~(mod((m2+2),3)),1);
%               
%               dat=datPre2(1:(m2(n2)+2));
%               
%           else      
%               
%         dat=datPre(1:(m(n)+2)); %% generate the sequence between start codon and stop codon  
%         end
%         
%        
%        else
%            continue;
%         end
%         
%        dat=datPre;
%         
%        N = floor(length(dat)/3);  %% amount of codons
%        
%         else
%         
%         dat=datPre;
%         N=length(dat)/3;
% 
%         ind = zeros (N, 3);
% 
%         ind = reshape (1:length(dat), 3, N)';
% 
%         pasteCodonPre = dat(ind);  %% generate codon sequence, portions of 3 nucleartides
% 
%         pasteCodonArray{i} = cellstr(pasteCodonPre); %% convert each unit type to cell for further manipulation
% 
%         end
%     end
% 
%     end


%%%%%%%%%%%%%%%%%%%%%%%only divide 3%%%%%%%%%%%%%%%%%%%%%

   
    fileID = fopen(filename);
    tt = textscan(fileID,'%s','delimiter',';');
    fclose(fileID);

    xx = tt{1};
    l = length(xx);
    datRaw = xx(2:2:l,:);   

    for i=1:length(datRaw)
        
        dat=datRaw{i};
        
        if mod(length(dat),3)==0
            
        N=length(dat)/3;

        ind = zeros (N, 3);

        ind = reshape (1:length(dat), 3, N)';

        pasteCodonPre = dat(ind);  %% generate codon sequence, portions of 3 nucleartides

        pasteCodonArray{i} = cellstr(pasteCodonPre); %% convert each unit type to cell for further manipulation

        else
        pasteCodonArray{i}=NaN;
        
        clear dat
        end
    end
end

    

