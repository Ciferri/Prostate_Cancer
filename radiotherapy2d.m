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

fid = fopen('tissueDim.dat', 'r');
dim = fscanf(fid, '%i');
nrow = dim(1);
ncol = dim(2);
nlayer = dim(3);
fclose(fid);

n3 = size(out, 1);

figure(1)

for i = 1:n3
  matrix(:,:,i) = reshape(out(i, :), ncol, nrow);
  img = uint16(matrix(:, :, i)');
  imshow(img, cmap)
  img = imresize(img, 4, 'nearest');
  path = ['img/figure' num2str(i, '%04d') '.png'];
  imwrite(img, cmap, path)
  pause(0.1)
end
