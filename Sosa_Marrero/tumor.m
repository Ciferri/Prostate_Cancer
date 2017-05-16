clear all
close all
load('out.dat');
cmap=[0.5 0.5 0.5; 
      1 1 1; 
      0 0 1;
      1 0 0;
      0 0 0
      ];
n1=94;
n2=57;
n3=size(out,1);
figure(1)
for i=1:n3
matrix(:,:,i)=reshape(out(i,:),n1,n2);
img=uint16(matrix(:,:,i)');
imshow(img,cmap)
path= ['/home/carlos/Documentos/Stage/v3/M2SLv0.1_EP/M2SLv0.1_EP/img/figure' num2str(i,'%04d') '.png'];
saveas(1,path);
pause(0.05)
end
