function diffusion()
  close all;
  
  x=107;
  y=66;
  
  palier1=0;
  palier2=0;
  palier3=0;
  palier4=0;
  palier5=0;
  palier6=0;
  palier7=0;
  palier8=0;
  palier9=0;
  palier10=0;
  
  load('out.dat');
  for i=1:50
    matrix(:,:,i)=reshape(out(i,:),x,y);
  end
    
  for i=1:50   
    figure(1);
    A = matrix(:,:,i);
    B = imrotate(A,-90);
    C = flip(B,2);
    imagesc(C);
    #colormap('gray');
    img = ['/home/ciferri/Bureau/stage/branch_github/Ciferri/img/impVes12/figure' num2str(i,'%04d') '.png'];
    saveas(1,img);
  end  
  
  for i=1:x  
    for j=1:y 
      A = matrix(i,j,50);
      if(A<1)
        palier1 += 1;
      elseif(A>=1 && A<5)
        palier2 +=1;
      elseif(A>=5 && A<10)
        palier3 +=1;
      elseif(A>=10 && A<15)
        palier4 +=1;
      elseif(A>=15 && A<20)
        palier5 +=1;
      elseif(A>=20 && A<25)
        palier6 +=1;
      elseif(A>=25 && A<30)
        palier7 +=1;
      elseif(A>=30 && A<35)
        palier8 +=1;
      elseif(A>=35 && A<40)
        palier9 +=1;
      else
        palier10 += 1;
        
      end
    end
  end
  inf1 = (palier1*58)/(x*y)
  sup1inf5 = (palier2*58)/(x*y)
  sup5inf10 = (palier3*58)/(x*y)
  sup10inf15 = (palier4*58)/(x*y)
  sup15inf20 = (palier5*58)/(x*y)
  sup20inf25 = (palier6*58)/(x*y)
  sup25inf30 = (palier7*58)/(x*y)
  sup30inf25 = (palier8*58)/(x*y)
  sup35inf40 = (palier9*58)/(x*y)
  palier10
  palier1
  palier2
  palier3
  palier4
  palier5
  palier6
  palier7
  palier8
  palier9
  palier10
  sizeTissue = x*y
endfunction
