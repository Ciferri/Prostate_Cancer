function impVes()
  img=imread("InitVes.png");
  fid=fopen('inVes.dat','w');
  for i=1:63
    for j=1:98
      if img(i,j) > 0
        fprintf(fid,'%f\n',1);
      else
        fprintf(fid,'%f\n',0);
      endif
    end
  end
  fclose(fid);
endfunction