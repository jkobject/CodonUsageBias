legendName={'(-2,0]','(-4,-2]','(-6,-4]','(-8,-6]','(-10,-8]','(-30,-10]','(-50,-30]','(-70,-50]','(-100,-70]'};

for i=1:9
figure
plot(N{1,i},'o');
xlabel('indices');
ylabel('subsequence length');
legend(legendName{i});
title('Arg, Sah6,subsequence length distribution corresponding to log(pr/pm) value shown as legend');
end

figure
plot(x16,y16,'o');
xlabel('midpoint of bin');
ylabel('log of frequency');
title('Arg, Sah6,(Pr/Pm) probability ratio distribution before replacement');
