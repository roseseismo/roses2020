      program likelihood
      implicit none
c
      real*8 hh,kk,pr,sori,var,tc_err,tc_exp,vv_err,vv_exp
c
      open(11,file='TA.M11A.xyz',status='old')
      rewind(11)
      open(12,file='TA.M11A.L.xyz',status='unknown')
      rewind(12)
c
c Expected values and one-sigma uncertainties from PyKrige:
      vv_exp=1.8070d0
      vv_err=0.03257d0
      tc_exp=33.2889d0
      tc_err=3.7293d0
c
c Read in the stack data:
1     read(11,*,end=2)kk,hh,sori
        var=((kk-vv_exp)/vv_err)**2+((hh-tc_exp)/tc_err)**2      ! Assume uncorrelated, normal
        pr=dexp(-var/2d0)
        write(12,*),kk,hh,pr,pr*sori
      goto1
2     close(11)
      close(12)
c
      end
