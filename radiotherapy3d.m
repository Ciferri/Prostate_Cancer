clear all
close all
pkg load image
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
for i=1:3:3*(n3-1)
matrix(:,:,i)=reshape(out(i,1:5358),n1,n2);
img=uint16(matrix(:,:,i)');
imshow(img,cmap)
img=imresize(img,4,'nearest');
path= ['img/figure' num2str(i,'%04d') '.png'];
imwrite(img,cmap,path)
pause(0.1)

matrix(:,:,i+1)=reshape(out(i,5359:10716),n1,n2);
img=uint16(matrix(:,:,i+1)');
imshow(img,cmap)
img=imresize(img,4,'nearest');
path= ['img/figure' num2str(i+1,'%04d') '.png'];
imwrite(img,cmap,path)
pause(0.1)

matrix(:,:,i+2)=reshape(out(i,10717:16074),n1,n2);
img=uint16(matrix(:,:,i+2)');
imshow(img,cmap)
img=imresize(img,4,'nearest');
path= ['img/figure' num2str(i+2,'%04d') '.png'];
imwrite(img,cmap,path)
pause(0.1)
end
