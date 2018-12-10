%% for synonymous codon of 2 and 3
% % AveEntropy2f801t1500=dlmread('AveEntropy2f801t1500.txt');
% % AveEntropy2(801:1500)=(AveEntropy2f801t1500)';
% % save 'AveEntropy2f.mat' AveEntropy2
% 
% AveEntropy3f601t800=dlmread('AveEntropy3f601t800.txt');
% AveEntropy3f801t1000=dlmread('AveEntropy3f801t1000.txt');
% 
% AveEntropy3(401:600)=(AveEntropy3f401t600)';
% AveEntropy3(601:800)=(AveEntropy3f601t800)';
% AveEntropy3(801:1000)=(AveEntropy3f801t1000)';
% 
% save 'AveEntropy3f.mat' AveEntropy3
% 
% % may not good enough
% samp1Mean6f201t499=zeros(1,299);
% for i=201:499
% filename1=['p1Sample6o',num2str(i),'.txt']; %% read from each txt file
% sampMt=dlmread(filename1);
% samp1Mean6f201t499(i-200)=mean(sampMt(1000:1000000)); %% samp1 of Mean6f201t499 formed here
% end
% 
% samp2Mean6f201t300=zeros(1,100);
% for j=201:300
% % m2=matfile('p22Sample6f201t300.mat');
% samp26f201t300=p22Sample6{j};
% samp2Mean6f201t300(j-200)=mean(samp26f201t300(1000:1000000));
% end
% 
% samp2Mean6f301t400=zeros(1,100);
% for u=301:400
% % m3m3=matfile('p33Sample6f301t400.mat');
% samp26f301t400=p33Sample6{u};
% samp2Mean6f301t400(u-300)=mean(samp26f301t400(1000:1000000));
% end
% 
% samp2Mean6f401t499=zeros(1,99);
% for v=401:499
% % m4=matfile('p44Sample6f401t500.mat');
% samp26f401t499=p44Sample6{v};
% samp2Mean6f401t499(v-400)=mean(samp26f401t499(1000:1000000));
% end
% 
% samp2Mean6f201t499(1:100)=samp2Mean6f201t300;
% samp2Mean6f201t499(101:200)=samp2Mean6f301t400;
% samp2Mean6f201t499(201:299)=samp2Mean6f401t499;  %% samp2 of Mean6f201t499 formed here
%% above for 6 synonymous may not good enough

samp1Mean6f201t300=zeros(1,100); %% this better method combine different part of the sample
for i=201:300
filename1=['p1Sample6o',num2str(i),'.txt'];
sampMt=dlmread(filename1);
samp362=p22Sample6{i};
sampT2(1:500000)=sampMt(100001:600000);  
sampT2(500001:1000000)=samp362(100001:600000);
samp1Mean6f201t300(i-200)=mean(sampT2);
clear sampMt sampT2 samp362
end


samp1Mean6f301t400=zeros(1,100);
for i=301:400
filename1=['p1Sample6o',num2str(i),'.txt'];
sampMt=dlmread(filename1);
samp363=p33Sample6{i};
sampT3(1:500000)=sampMt(100001:600000);  
sampT3(500001:1000000)=samp363(100001:600000);
samp1Mean6f301t400(i-300)=mean(sampT3);
clear sampMt sampT3 samp363
end

samp1Mean6f401t499=zeros(1,99);
for i=401:499
filename1=['p1Sample6o',num2str(i),'.txt'];
sampMt=dlmread(filename1);
samp364=p44Sample6{i};
sampT4(1:500000)=sampMt(200001:700000);  
sampT4(500001:1000000)=samp364(200001:700000);
samp1Mean6f401t499(i-400)=mean(sampT4);
clear sampMt sampT4 samp364
end

samp1Mean6f201t499(1:100)=samp1Mean6f201t300;
samp1Mean6f201t499(101:200)=samp1Mean6f301t400;
samp1Mean6f201t499(201:299)=samp1Mean6f401t499;   %% sample3 of mean6 from 201 to 499 formed here


samp2Mean6f201t300=zeros(1,100);
for i=201:300
filename1=['p1Sample6o',num2str(i),'.txt'];
sampMt=dlmread(filename1);
samp362=p22Sample6{i};
sampT2(1:500000)=sampMt(300001:800000);  
sampT2(500001:1000000)=samp362(300001:800000);
samp2Mean6f201t300(i-200)=mean(sampT2);
clear sampMt sampT2 samp362
end


samp2Mean6f301t400=zeros(1,100);
for i=301:400
filename1=['p1Sample6o',num2str(i),'.txt'];
sampMt=dlmread(filename1);
samp363=p33Sample6{i};
sampT3(1:500000)=sampMt(300001:800000);  
sampT3(500001:1000000)=samp363(300001:800000);
samp2Mean6f301t400(i-300)=mean(sampT3);
clear sampMt sampT3 samp363
end

samp2Mean6f401t499=zeros(1,99);
for i=401:499
filename1=['p1Sample6o',num2str(i),'.txt'];
sampMt=dlmread(filename1);
samp364=p44Sample6{i};
sampT4(1:500000)=sampMt(300001:800000);  
sampT4(500001:1000000)=samp364(300001:800000);
samp2Mean6f401t499(i-400)=mean(sampT4);
clear sampMt sampT4 samp364
end

samp2Mean6f201t499(1:100)=samp2Mean6f201t300;
samp2Mean6f201t499(101:200)=samp2Mean6f301t400;
samp2Mean6f201t499(201:299)=samp2Mean6f401t499;   %% sample3 of mean6 from 201 to 499 formed here


samp1Mean6f101t500(101:200)=samp1Mean6f101t200;
samp1Mean6f101t500(201:499)=samp1Mean6f201t499;
samp1Mean6f101t500(500)=mean(p1Sample6(1000000:2000000));

samp2Mean6f101t500(101:200)=samp2Mean6f101t200;
samp2Mean6f101t500(201:499)=samp2Mean6f201t499;
samp2Mean6f101t500(500)=mean(p2Sample6(3000000:4000000));

samp3Mean6f201t300=zeros(1,100);
for i=201:300
filename1=['p1Sample6o',num2str(i),'.txt'];
sampMt=dlmread(filename1);
samp362=p22Sample6{i};
sampT2(1:500000)=sampMt(200001:700000);  
sampT2(500001:1000000)=samp362(200001:700000);
samp3Mean6f201t300(i-200)=mean(sampT2);
clear sampMt sampT2 samp362
end


samp3Mean6f301t400=zeros(1,100);
for i=301:400
filename1=['p1Sample6o',num2str(i),'.txt'];
sampMt=dlmread(filename1);
samp363=p33Sample6{i};
sampT3(1:500000)=sampMt(200001:700000);  
sampT3(500001:1000000)=samp363(200001:700000);
samp3Mean6f301t400(i-300)=mean(sampT3);
clear sampMt sampT3 samp363
end

samp3Mean6f401t499=zeros(1,99);
for i=401:499
filename1=['p1Sample6o',num2str(i),'.txt'];
sampMt=dlmread(filename1);
samp364=p44Sample6{i};
sampT4(1:500000)=sampMt(200001:700000);  
sampT4(500001:1000000)=samp364(200001:700000);
samp3Mean6f401t499(i-400)=mean(sampT4);
clear sampMt sampT4 samp364
end

samp3Mean6f201t499(1:100)=samp3Mean6f201t300;
samp3Mean6f201t499(101:200)=samp3Mean6f301t400;
samp3Mean6f201t499(201:299)=samp3Mean6f401t499;   %% sample3 of mean6 from 201 to 499 formed here

%%calculate mean and standard error
sampMeanWhole=[samp1Mean6f101t500;samp2Mean6f101t500;samp3Mean6f101t500];
% sampMeanWhole=[samp1Mean6f101t500;samp2Mean6f101t500;samp3Mean6f101t500;samp4Mean6f101t500;samp5Mean6f101t500];
MeanOSmean=mean(sampMeanWhole);
StandEr=std(sampMeanWhole);

% TblSampleMeanAnalysis = table((samp1Mean6f101t500)',(samp2Mean6f101t500)',(samp3Mean6f101t500)',...
%     (samp4Mean6f101t500)',(samp5Mean6f101t500)',(MeanOSmean)',(StandEr)',...
%     'VariableName',{'mean1 of sample1' 'mean2 of sample2' 'mean3 of sample3' ...
%     'mean4 of sample4' 'mean5 of sample5' 'mean of means' 'standard deviation of means'});

%%generate table and write to txt file
TblSampleMeanAnalysis = table((samp1Mean6f101t500)',(samp2Mean6f101t500)',(samp3Mean6f101t500)',...
    (MeanOSmean)',(StandEr)',...
     'VariableName',{'mean1Sample1' 'mean2Sample2' 'meanSsample3' ...
     'meanOFmeans' 'StdOfMeans'});

writetable(TblSampleMeanAnalysis)

%%% process the sample value log(pi) into entropy cost log(pi/pmax)
% for i=101:500
% pmax=Efor(6,i);
% p=exp(MeanOSmean(i));
% temp(i-100)=log(p/pmax);
% end


