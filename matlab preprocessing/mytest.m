function [] =mytest(l)

fid=fopen('mytest.txt','w');
fprintf(fid,'%d,',l);
fclose(fid);

end