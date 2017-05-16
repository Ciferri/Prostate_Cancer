img=imread("InitVes.png");
fid=fopen('inVes.dat','w');
for i=1:57
for j=1:94
fprintf(fid,'%f\n',img(i,j));
end
end
fclose(fid);