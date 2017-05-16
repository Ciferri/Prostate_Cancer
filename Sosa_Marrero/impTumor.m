img=imread("InitTumor.png");
fid=fopen('inTumor.dat','w');
for i=1:57
for j=1:94
fprintf(fid,'%f\n',img(i,j));
end
end
fclose(fid);