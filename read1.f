	double precision lambda,re1,re2,re3,im1,im2,im3

        open(1,file="ShettleSandDust.dat")

        open(2,file="Sand.dat")
        open(3,file="Dust.dat")

        do 100 i=1,68
          read(1,*)lambda,re1,im1,re2,im2,re3,im3
          write(2,*)lambda,(re1+re2)/2d0,(im1+im2)/2d0
          write(3,*)lambda,re3,im3
100    continue
       stop
       end
