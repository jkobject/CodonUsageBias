setSynonymousCodonTable;

% % GeneInterestList=getHomologyList('saccharomyces_arboricola_h_6Gene.fa');%%%%generate homoList which contain interested genes
% % indChosen=randi(3659,[500,1]);
% % homoListSel=GeneInterestList(indChosen,:);

% % % homoList={'YGL103W','YOL127W','YLR075W','YBL092W','YOR063W','YLR340W','YOR369C',...
% % %     'YDR064W','YOL040C','YGL123W','YPL267W','YIL158W','YCL029C','YNL166C','YOR026W',...
% % %     'YCL014W','YLR175W','YCR002C','YLR215C','YKL022C','YCL040W','YFR053C','YGL059W',...
% % %     'YGL253W','YGR192C','YHR044C','YJL155C','YKL127W','YMR278W','YNL241C','YJL200C',...
% % %     'YGR204W','YAL012W','YIL116W','YMR108W','YLR451W','YDL182W','YNL277W','YDR007W',...
% % %     'YOR323C','YFR028C','YLR272C','YBR135W','YMR199W','YAL040C','YJL164C','YGL240W'...
% % %     'YML027W','YDR201W','YDR293C'};   %%homology sourse


% fileID=fopen('homoList5818Rand.txt');
% tt=textscan(fileID,'%s','Delimiter',',');
% homoList=tt{1,1};
% fclose(fileID);

homoList={'YBR221C'};

% for homoc=1:length(homoList)

for homoc=1:1 %%homology counter
    
    stringName=homoList{homoc};
    
    [homoSpecies,homoGene,l] = getHomoInfor(stringName);%%%call function 'getHomoInfor' to retrieve homology information
    
    if l~=0
    
    fileID1=fopen([stringName,'homologyGlu.txt'],'a');
    fprintf(fileID1,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID2=fopen([stringName,'homologyHis.txt'],'a');
    fprintf(fileID2,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID3=fopen([stringName,'homologyGln.txt'],'a');
    fprintf(fileID3,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID4=fopen([stringName,'homologyLys.txt'],'a');
    fprintf(fileID4,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID5=fopen([stringName,'homologyAsp.txt'],'a');
    fprintf(fileID5,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID6=fopen([stringName,'homologyPhe.txt'],'a');
    fprintf(fileID6,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID7=fopen([stringName,'homologyTyr.txt'],'a');
    fprintf(fileID7,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID8=fopen([stringName,'homologyCys.txt'],'a');
    fprintf(fileID8,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID9=fopen([stringName,'homologyAsn.txt'],'a');
    fprintf(fileID9,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID10=fopen([stringName,'homologyIle.txt'],'a');
    fprintf(fileID10,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID11=fopen([stringName,'homologyGly.txt'],'a');
    fprintf(fileID11,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID12=fopen([stringName,'homologyPro.txt'],'a');
    fprintf(fileID12,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID13=fopen([stringName,'homologyVal.txt'],'a');
    fprintf(fileID13,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID14=fopen([stringName,'homologyThr.txt'],'a');
    fprintf(fileID14,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID15=fopen([stringName,'homologyAla.txt'],'a');
    fprintf(fileID15,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID16=fopen([stringName,'homologyArg.txt'],'a');
    fprintf(fileID16,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID17=fopen([stringName,'homologySer.txt'],'a');
    fprintf(fileID17,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID18=fopen([stringName,'homologyLeu.txt'],'a');
    fprintf(fileID18,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID19=fopen([stringName,'homologyError.txt'],'a');
    
    
        for j=1:l
            try
                speciesName=homoSpecies{j};
                fileName1=sprintf('%s',[speciesName,'.fa']);%% find corresponding species file
                fileName2=sprintf('%s',[speciesName,'Gene.fa']);%% only gene names in such file
                geneNamesWhole=getHomologyList(fileName2); %% list all the gene names in the interesed species
                geneName=homoGene{j};
                
                id=find(contains(geneNamesWhole,geneName));
                pasteCodonWhole=getCodonSequence(fileName1); %%%%% find interested homology
                pasteCodon=pasteCodonWhole{id};
                
% % %                 %%%%Glu%%%%%%
% % %                 [NNe,~,Ep]=GluAminoAcidH(pasteCodon);
% % %                 E=-log(Ep)/NNe;
% % %                 if NNe<=800
% % %                     Plocation = entropy4Location(2,NNe,E);
% % %                 else
% % %                     Plocation=NaN;
% % %                 end
% % %                 
% % %                 fprintf(fileID1,'%s,%s,%f,%f,%i\n',speciesName,geneName,E,Plocation,NNe);
% % %                 
% % %           %%%%%%%%%%%%%%%His%%%%%%
%                 [NNh,~,Hp]=HisAminoAcidH(pasteCodon);
%                 H=-log(Hp)/NNh;
%                 if NNh<=800
%                     Plocation = entropy4Location(2,NNh,H);
%                 else
%                     Plocation=NaN;
%                 end
%                 
%                 fprintf(fileID2,'%s,%s,%f,%f,%i\n',speciesName,geneName,H,Plocation,NNh);
%                 
% % %                 %%%%%Gln%%%%%%
% % %                 [NNq,~,Qp]=GlnAminoAcidH(pasteCodon);
% % %                 Q=-log(Qp)/NNq;
% % %                 if NNq<=800
% % %                     Plocation = entropy4Location(2,NNq,Q);
% % %                 else
% % %                     Plocation=NaN;
% % %                 end
% % %                 
% % %                 fprintf(fileID3,'%s,%s,%f,%f,%i\n',speciesName,geneName,Q,Plocation,NNq);
% % %                 
% % %                 %%%%Lys%%%%%%
% % %                 [NNk,~,Kp]=LysAminoAcidH(pasteCodon);
% % %                 K=-log(Kp)/NNk;
% % %                 if NNk<=800
% % %                     Plocation = entropy4Location(2,NNk,K);
% % %                 else
% % %                     Plocation=NaN;
% % %                 end
% % %                 
% % %                 fprintf(fileID4,'%s,%s,%f,%f,%i\n',speciesName,geneName,K,Plocation,NNk);
% % %                 
% % %                 %%%%%Asp%%%%%%
% % %                 [NNd,~,Dp]=AspAminoAcidH(pasteCodon);
% % %                 D=-log(Dp)/NNd;
% % %                 if NNd<=800
% % %                     Plocation = entropy4Location(2,NNd,D);
% % %                 else
% % %                     Plocation=NaN;
% % %                 end
% % %                 
% % %                 fprintf(fileID5,'%s,%s,%f,%f,%i\n',speciesName,geneName,D,Plocation,NNd);
% % %                 
% % %                 %%%%%Phe%%%%%%
% % %                 [NNf,~,Fp]=PheAminoAcidH(pasteCodon);
% % %                 F=-log(Fp)/NNf;
% % %                 if NNf<=800
% % %                     Plocation = entropy4Location(2,NNf,F);
% % %                 else
% % %                     Plocation=NaN;
% % %                 end
% % %                 
% % %                 fprintf(fileID6,'%s,%s,%f,%f,%i\n',speciesName,geneName,F,Plocation,NNf);
% % %                 
% % %                 %%%%%Tyr%%%%%%
% % %                 [NNy,~,Yp]=TyrAminoAcidH(pasteCodon);
% % %                 Y=-log(Yp)/NNy;
% % %                 if NNy<=800
% % %                     Plocation = entropy4Location(2,NNy,Y);
% % %                 else
% % %                     Plocation=NaN;
% % %                 end
% % %                 
% % %                 fprintf(fileID7,'%s,%s,%f,%f,%i\n',speciesName,geneName,Y,Plocation,NNy);
% % %                 
                %%%%%Cys%%%%%%
                [NNc,~,Cp]=CysAminoAcidH(pasteCodon);
                C=-log(Cp)/NNc;
                if NNc<=800
                    Plocation = entropy4Location(2,NNc,C);
                else
                    Plocation=NaN;
                end
                
                fprintf(fileID8,'%s,%s,%f,%f,%i\n',speciesName,geneName,C,Plocation,NNc);
                
                %%%%%%%%%Asn%%%%%%
%                 [NNn,~,Np]=AsnAminoAcidH(pasteCodon);
%                 N=-log(Np)/NNn;
%                 if NNn<=800
%                     Plocation = entropy4Location(2,NNn,N);
%                 else
%                     Plocation=NaN;
%                 end
%                 
%                 fprintf(fileID9,'%s,%s,%f,%f,%i\n',speciesName,geneName,N,Plocation,NNn);
%                 
                %%%%%%%%%%Ile%%%%%%
%                 [NNi,~,Ip]=IleAminoAcidH(pasteCodon);
%                 I=-log(Ip)/NNi;
%                 if NNi<=400
%                     Plocation = entropy4Location(3,NNi,I);
%                 else
%                     Plocation=NaN;
%                 end
%                 
%                 fprintf(fileID10,'%s,%s,%f,%f,%i\n',speciesName,geneName,I,Plocation,NNi);
%                 
%                 %%%%%%%%%%Gly%%%%%%%
%                 [NNg,~,Gp]=GlyAminoAcidH(pasteCodon);
%                 G=-log(Gp)/NNg;
%                 if NNg<=300
%                     Plocation = entropy4Location(4,NNg,G);
%                 else
%                     Plocation=NaN;
%                 end
%                 
%                 fprintf(fileID11,'%s,%s,%f,%f,%i\n',speciesName,geneName,G,Plocation,NNg);
% %                 
% % %                 %%%%%%%%%%Pro%%%%%%%
%                 [NNp,~,Pp]=ProAminoAcidH(pasteCodon);
%                 P=-log(Pp)/NNp;
%                 if NNp<=300
%                     Plocation = entropy4Location(4,NNp,P);
%                 else
%                     Plocation=NaN;
%                 end
%                 
%                 fprintf(fileID12,'%s,%s,%f,%f,%i\n',speciesName,geneName,P,Plocation,NNp);
% % %                 
% % %                 %%%%%%%%Val%%%%%%%
% % %                 [NNv,~,Vp]=ValAminoAcidH(pasteCodon);
% % %                 V=-log(Vp)/NNv;
% % %                 if NNv<=300
% % %                     Plocation = entropy4Location(4,NNv,V);
% % %                 else
% % %                     Plocation=NaN;
% % %                 end
% % %                 
% % %                 fprintf(fileID13,'%s,%s,%f,%f,%i\n',speciesName,geneName,V,Plocation,NNv);
% % %                 
% % %                 %%%%%%%%%Thr%%%%%%%
% % %                 [NNt,~,Tp]=ThrAminoAcidH(pasteCodon);
% % %                 T=-log(Tp)/NNt;
% % %                 if NNt<=300
% % %                     Plocation = entropy4Location(4,NNt,T);
% % %                 else
% % %                     Plocation=NaN;
% % %                 end
% % %                 
% % %                 fprintf(fileID14,'%s,%s,%f,%f,%i\n',speciesName,geneName,T,Plocation,NNt);
% % %                 
% % %           %%%%%%%%%Ala%%%%%%%
% % %                 [NNa,~,Ap]=AlaAminoAcidH(pasteCodon);
% % %                 A=-log(Ap)/NNa;
% % %                 if NNa<=300
% % %                     Plocation = entropy4Location(4,NNa,A);
% % %                 else
% % %                     Plocation=NaN;
% % %                 end
% % %                 
% % %                 fprintf(fileID15,'%s,%s,%f,%f,%i\n',speciesName,geneName,A,Plocation,NNa);
% % %                 
                %%%%%%%%%%%%%%Arg%%%%%%%
%                 [NNr,~,Rp]=ArgAminoAcidH(pasteCodon);
%                 R=-log(Rp)/NNr;
%                 if NNr<=400
%                     Plocation = entropy4Location(6,NNr,R);
%                 else
%                     Plocation=NaN;
%                 end
%                 
%                 fprintf(fileID16,'%s,%s,%f,%f,%i\n',speciesName,geneName,R,Plocation,NNr);
%                 
% % % % %                 %%%%%%%%Ser%%%%%%%
% % % %                 [NNs,~,Sp]=SerAminoAcidH(pasteCodon);
% % % %                 S=-log(Sp)/NNs;
% % % %                 if NNs<=400
% % % %                     Plocation = entropy4Location(6,NNs,S);
% % % %                 else
% % % %                     Plocation=NaN;
% % % %                 end
% % % %                 
% % % %                 fprintf(fileID17,'%s,%s,%f,%f,%i\n',speciesName,geneName,S,Plocation,NNs);
% % % %                 
% % % %                 %%%%Leu%%%%%%%
% % % %                 [NNl,~,Lp]=LeuAminoAcidH(pasteCodon);
% % % %                 L=-log(Lp)/NNl;
% % % %                 if NNl<=400
% % % %                     Plocation = entropy4Location(6,NNl,L);
% % % %                 else
% % % %                     Plocation=NaN;
% % % %                 end
% % % %                 
% % % %                 fprintf(fileID18,'%s,%s,%f,%f,%i\n',speciesName,geneName,L,Plocation,NNl);
                
            catch
                fprintf(fileID19,'%s,%s,%s,%s\n',[speciesName,":",geneName,':error!']);
            end
            
            clear pasteCodon
        end
        
    fclose(fileID1);
    fclose(fileID2);
    fclose(fileID3);
    fclose(fileID4);
    fclose(fileID5);
    fclose(fileID6);
    fclose(fileID7);
    fclose(fileID8);
    fclose(fileID9);
    fclose(fileID10);
    fclose(fileID11);
    fclose(fileID12);
    fclose(fileID13);
    fclose(fileID14);
    fclose(fileID15);
    fclose(fileID16);
    fclose(fileID17);
    fclose(fileID18);
    fclose(fileID19);
    end
end


