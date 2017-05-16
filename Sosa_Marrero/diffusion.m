function diffusion()
  close all;
  load('out.dat');
  
  for i=1:10
    matrix(:,:,i)=reshape(out(i,:),5,5);
  end
    
  for i=1:10,   
    figure(1);
    imagesc(matrix(:,:,i));
    %toto=matrix(:,:,i);
    colormap('gray');
    img = ['/home/carlos/Documentos/Stage/v3/M2SLv0.1_EP/M2SLv0.1_EP/img/figure' num2str(i,'%04d') '.png'];
    saveas(1,img);
    %imwrite(toto,img);
  end
endfunction
