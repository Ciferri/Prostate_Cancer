clear all
close all
img=imread("initTissues/Tissu P14 01685 7-5-cd31 97x57.png");
nrow=57;
ncol=97;
img=imresize(img,[nrow,ncol]);
imgB=img(:,:,3);
initVes=imgB>=0 & imgB<80;
fid=fopen('inVes.dat','w');
count=0;
for i=1:nrow
  for j=1:ncol
    fprintf(fid,'%f\n',initVes(i,j));
    if(initVes(i,j)==1)
    count++;
    endif
  end
end
count
fclose(fid);
initTum=imgB>80 &imgB<200;
fid=fopen('inTum.dat','w');
count=0;
for i=1:nrow
  for j=1:ncol
    fprintf(fid,'%f\n',initTum(i,j));
    if(initTum(i,j)==1)
    count++;
    endif
  end
end
count
fclose(fid);

initPO2=0.77*ones(nrow,ncol);
fid=fopen('inPO2.dat','w');
for i=1:nrow
  for j=1:ncol
      fprintf(fid,'%f\n',initPO2(i,j));
  end
end
fclose(fid);