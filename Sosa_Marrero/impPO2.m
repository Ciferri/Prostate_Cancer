img=3.5*ones(57,94);
fid=fopen('inPO2.dat','w');
for i=1:57
for j=1:94
fprintf(fid,'%f\n',img(i,j));
end
end
fclose(fid);