img=zeros(57,94);
img(29,47)=1;
fid=fopen('inTumor2.dat','w');
for i=1:57
for j=1:94
fprintf(fid,'%f\n',img(i,j));
end
end
fclose(fid);