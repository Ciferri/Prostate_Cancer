clear all
close all
img=imread("initTissues/Tissu P14 26671 4-1-cd31 94x58.png");
nrow=58;
ncol=94;

img=imresize(img,[nrow,ncol]);
imgB=img(:,:,3);

initVes=imgB>=0 & imgB<80;
fid=fopen('inVes.dat','w');
for i=1:nrow
  for j=1:ncol
    fprintf(fid,'%f\n',initVes(i,j));
  end
end
fclose(fid);

initTum=imgB>80 &imgB<200;
fid=fopen('inTum.dat','w');
for i=1:nrow
  for j=1:ncol
    fprintf(fid,'%f\n',initTum(i,j));
  end
end
fclose(fid);

load('PO2.dat');
n3=size(PO2,1);
fid=fopen('inPO2.dat','w');
initPO2=reshape(PO2(n3,:),ncol,nrow)';
for i=1:nrow
  for j=1:ncol
      fprintf(fid,'%f\n',initPO2(i,j));
  end
end
fclose(fid);