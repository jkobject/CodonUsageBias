types={'ribosomal protein','transcription factor','lucose metabolic process',...
    'Amino acid biosynthetic process','Cell Cycle'};

homoList={'YGL103W','YOL127W','YLR075W','YBL092W','YOR063W','YLR340W','YOR369C','YDR064W','YOL040C','YGL123W',...
    'YPL267W','YIL158W','YCL029C','YNL166C','YOR026W','YCL014W','YLR175W','YCR002C','YLR215C','YKL022C',...
    'YCL040W','YFR053C','YGL059W','YGL253W','YGR192C','YHR044C','YJL155C','YKL127W','YMR278W','YNL241C',...
    'YJL200C','YGR204W','YAL012W','YIL116W','YMR108W','YLR451W','YDL182W','YNL277W','YDR007W','YOR323C',...
    'YFR028C','YLR272C','YBR135W','YMR199W','YAL040C','YJL164C','YGL240W','YML027W','YDR201W','YDR293C'};

% AAname={'Glu','His','Gln','Phe','Tyr','Cys','Asn','Lys','Asp','Ile',...
%     'Pro','Thr','Ala','Val','Gly','Arg','Ser','Leu'};
AAname={'Asp','Ile','Gly','Arg'};


%%%%%per plot: per amino acid/ all groups

for an=1:length(AAname)  %%%amino acid counter: an
    aaName=AAname{an};
        [XfPlot,grp]=getHomoPlot3(homoList,aaName);
        figure
        boxplot(XfPlot,grp,'LabelOrientation','inline');
        txt = findobj(gca,'Type','text');   %%%% set xLabel font size
        set(txt(1:end),'VerticalAlignment', 'Middle','FontSize',8);
        title(['50 Artificial homology analysis: Amino Acid: ',aaName]);
end

%%%%per plot: per amino acid/per group 
% for an=1:length(AAname)  %%%amino acid counter: an
%     aaName=AAname{an};
% 
%     for tt=0:4  %%group counter:tt
%         [XfPlot,grp]=getHomoPlot(tt,homoList,aaName);
%         figure
%         boxplot(XfPlot,grp);
%         title([types{tt+1},': Amino Acid: ',aaName]);
%     end
% 
% end


%%%%per plot: 18 amino acids/ per group
% for tt=0:4
%     [XfPlot,grp]=getHomoPlot2(tt,homoList,AAname);
%     
%     figure
%     boxplot(XfPlot,grp);
%     title(types{tt+1});
% end
