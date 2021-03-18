***********************************************     
*     delta term at nlo
***     
      double precision function delta1(lfr)
      implicit none

      double precision pi,cf,z2
      double precision lfr,lfac

c     lfr here is actually -lfh
      lfac = -lfr

      cf = 4d0 / 3d0
      pi = 3.1415926535897932385d0
      z2 = 1.6449340668482264365d0
      
      delta1 = 8.0d0 * z2 - 16.0d0 + 6.0d0 * lfac
      delta1 = cf * delta1/4.d0/pi

      return
      end

***     
*     subtraction term associated with + prescription at nlo
***     
      double precision function dtsub1(lfh,tauh)
      implicit none
      double precision cf,lfac,pi
      double precision ltmin,lfh,tauh

c     lfr here is actually -lfh
      lfac = -lfh
      pi = 3.1415926535897932385d0
      cf = 4.d0/3d0
      
      ltmin = dlog(1.0d0 - tauh)
      dtsub1 = 8.0d0 * ltmin**2  + 8.0d0*ltmin * lfac
      dtsub1 = cf * dtsub1/4.d0/pi

      return
      end

***
*     plus at nlo!
***
      double precision function dterms1(xt,lfh)
      implicit none

      double precision ca,cf,nf,pi
      double precision xt,lfh,lfac
      double precision dd0,dd1

c     lfh here is actually -lfh
      lfac = -lfh
      pi = 3.1415926535897932385d0
      
      ca = 3.d0
      cf = 4d0/3d0
      nf = 5d0
      dd0 = 1.d0/(1.0d0-xt)
      dd1 = dd0 * dlog(1.0d0-xt)
      dterms1 = 16.0d0 * dd1
      dterms1 = dterms1 + 8.0d0 * lfac * dd0
      dterms1 = cf * dterms1/4.d0/pi

      return
      end
      
***
*     regular terms at nlo in qq
***
      double precision function sbbq1(xt,lfh)
      implicit none

      double precision ca,cf,nf
      double precision xt,lfh,lfac
      double precision dd0,dd1,pi

c     lfr here is actually -lfh
      lfac = -lfh
      pi = 3.1415926535897932385d0
 
      ca = 3.d0
      cf = 4d0/3d0
      nf = 5d0
      
      sbbq1 = -8.0d0*( 1.0d0 + xt )*dlog(1.0d0-xt)
     .     -4.0d0*( 1.0d0 + xt**2 )*dlog(xt)/(1.0d0-xt)
      sbbq1 = sbbq1 - 4.0d0*( 1.0d0 + xt )*lfac
      sbbq1 = xt* cf * sbbq1/4.d0/pi

      return
      end
***
*     term in front of qg function at nlo
***
      double precision function sbg1(xt,lfh)
      implicit none

      double precision ca,tf
      double precision xt,lfh,lfac
      double precision z,zmin,log1mz,pi

c     lfh here is actually -lfh
      lfac = -lfh
      pi = 3.1415926535897932385d0
 
      tf = 0.500000000000000000000d0
      z      = xt
      zmin   = 1.0d0-z
      log1mz = dlog(zmin)
      sbg1   = 2.0d0*( z**2 + zmin**2 )*( 2.0d0*log1mz - dlog(z) ) +
     .     1.d0 + 6.0d0*z - 7.0d0*z**2
      sbg1 = sbg1 + 2.d0*( z**2 + zmin**2 )*lfac
      sbg1 = xt * tf * sbg1/4.d0/pi
      return
      end


c.......................................................................
c............  NNLO .................................................
c.......................................................................

***
*     Regular term in bb~
**


      double precision function deltabbqa(xx,lfh,lfr)
      implicit none
      
      double precision xx,lfh,lfr,z,lfac,lfren
      double precision CQQIV,CQQ,CQQIA,CQBBB,CQBACD
      double precision CQBBCD,CQBAB,CQBCAR,CQBNFR,CQBCFR
      double precision cax,cv,cd,sw,pi
      external CQQIV,CQQ,CQQIA,CQBBB,CQBACD
      external CQBBCD,CQBAB,CQBCAR,CQBNFR,CQBCFR

      sw  = 0.2228972225239183d0
      cax = 1.d0
      cv  = -1.d0 + 4.d0/3.d0*sw
      cd  = 1.d0 + cv**2
      pi = 3.1415926535897932385d0
 
c     lfh here is actually -lfh,while lfr = -lfh+lfr
      lfren = -lfh+lfr
      lfac = -lfh
      
      deltabbqa = xx*(cd*(CQBBB(xx,lfac,lfren) + 2.d0*CQQ(xx,lfac,lfren)
     ?     + 2.d0*CQBBCD(xx,lfac,lfren) )
     ?     - 2.d0*(cv**2)*CQQIV(xx,lfac,lfren)
     ?     + 2.d0*(cax**2)*CQQIA(xx,lfac,lfren)
     ?     + CQBAB(xx,lfac,lfren)
     ?     + cd*(CQBCFR(xx,lfac,lfren) + CQBCAR(xx,lfac,lfren)
     ?     + CQBNFR(xx,lfac,lfren) + 2.d0*CQBACD(xx,lfac,lfren)))

      deltabbqa = deltabbqa/cd
      deltabbqa = deltabbqa/(4.d0*pi)**2
      return
      end

c.......................................................................

      double precision function deltabbqf(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr
      deltabbqf = 0d0
      return
      end
***
*     Regular term in bg
**
      double precision function deltabga(xx,lfh,lfr)
      implicit none
       
      double precision xx,lfh,lfr,lfac,lfren
      double precision xm1,xp1,dlnx,dlxm1,dlxp1
      double precision dli2a,dli3a,zcom,s1mz,dimz,trimz
      double precision trimco,tripco
      double precision z2,z3,z4
      double precision ca, cf, nf,pi
      double precision partCA, partCF,scaleCA,scaleCF,scaleNF
      double precision dilog,trilog,wgplg
      external dilog,trilog,wgplg
c     lfh here is actually -lfh,while lfr = -lfh+lfr
      lfren = -lfh+lfr
      lfac = -lfh      
      pi = 3.1415926535897932385d0

      ca = 3.d0
      cf = 4.d0/3.d0
      nf = 5.d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)
      xp1   = 1.0d0 + xx
      dlxp1 = dlog(xp1)
      dli2a = dilog(xm1)
      dli3a = trilog(xm1)
      dimz  = wgplg(1,1,-xx)
      trimz = wgplg(2,1,-xx)
      zcom  = xm1/xp1
      tripco= wgplg(2,1,zcom)
      trimco= wgplg(2,1,-zcom)
      s1mz  = wgplg(1,2,xm1)

      partCA =
     .   dlnx**3 * (  - 3.d0 - 20.d0/3.d0*xx )
     . + dlnx**2*dlxm1 * ( 6.d0 - 4.d0*xx**2 + 28.d0*xx )
     . + dlnx**2*dlxp1 * (  - 6 - 12*xx**2 - 12*xx )
     . + dlnx**2 * ( 5.d0/2.d0 + 173.d0/3.d0*xx**2 + 2.0d0*xx )
     . + dlnx*dlxm1**2 * (  - 2.d0 + 12.d0*xx**2 - 44.d0*xx )
     . + dlnx*dlxm1*dlxp1 * ( 8.d0 + 16.d0*xx**2 + 16.d0*xx )
     . + dlnx*dlxm1 * (  - 10.d0 - 130.d0*xx**2 + 16.d0*xx )
     . + dlnx*dlxp1 * (  - 4.d0 - 16.d0*xx**2 - 20.d0*xx )
     . + dlnx*dli2a * ( 8.d0*xx**2 - 28.d0*xx )
     . + dlnx*dimz * (  - 8.d0 - 16.d0*xx**2 - 16.d0*xx )
     . + dlnx*z2 * (  - 16.d0*xx**2 + 40.d0*xx )
     . + dlnx * (  - 118.d0/3.d0 + 457.d0/9.d0*xx**2 + 16.d0/3.d0*xx )
     . + dlxm1**3 * (  - 13.d0/3.d0 - 26.d0/3.d0*xx**2 + 26.d0/3.d0*xx )
     . + dlxm1**2 * (  - 4.d0 + 154.d0/3.d0*xx**2 - 42.d0*xx
     .     - 16.d0/3.d0/xx )
     . + dlxm1*dli2a * (  - 28.d0 - 20.d0*xx**2 - 40.d0*xx )
     . + dlxm1*dimz * ( 8.d0 + 16.d0*xx**2 + 16.d0*xx )
     . + dlxm1*z2 * ( 16.d0 + 32.d0*xx**2 - 16.d0*xx )
      partCA = partCA
     . + dlxm1 * ( 70.d0/3.d0 - 74.d0/9.d0*xx**2
     .     - 25.d0/3.d0*xx - 88.d0/9.d0/xx )
     .     + dli2a * (  - 22.d0 - 88.d0/3.d0*xx**2 - 64.d0*xx
     .     - 32.d0/3.d0/xx )
     . + dimz * (  - 4.d0 - 16.d0*xx**2 - 20.d0*xx )
     . + z2 * (  - 10.d0 - 214.d0/3.d0*xx**2 + 56.d0*xx
     .     + 16.d0/3.d0/xx )
     . + dli3a * ( 30.d0 + 24.d0*xx**2 + 68.d0*xx )
     . + trimz * ( 4 + 8*xx**2 + 8*xx )
     . + trimco * ( 8.d0 + 16.d0*xx**2 + 16.d0*xx )
     . + tripco * (  - 8.d0 - 16.d0*xx**2 - 16.d0*xx )
     . + s1mz * (  - 36.d0 - 16.d0*xx**2 - 64.d0*xx )
     . + z3 * ( 2.d0 + 4.d0*xx**2 + 8.d0*xx )
     . - 539.d0/18.d0 - 1837.d0/54.d0*xx**2 + 613.d0/9.d0*xx
     . - 58.d0/27.d0/xx

      partCA  = - ca * partCA

* -- CF 
      partCF =
     .   dlnx**3 * ( 17.d0/6.d0 + 26.d0/3.d0*xx**2 - 17.d0/3.d0*xx )
     . + dlnx**2*dlxm1 * (  - 12.d0 - 40.d0*xx**2 + 24.d0*xx )
     . + dlnx**2 * ( 11.d0/4.d0 + xx**2 - 15.d0*xx )
     . + dlnx*dlxm1**2 * ( 21.d0 + 66.d0*xx**2 - 42.d0*xx )
     . + dlnx*dlxm1 * (  - 14.d0 - 96.d0*xx**2 + 96.d0*xx )
     . + dlnx*dlxp1 * (  - 8.d0 + 24.d0*xx**2 + 16.d0*xx )
     . + dlnx*dli2a * (  - 2.d0 + 4.d0*xx )
     . + dlnx*dimz * ( 8.d0 + 16.d0*xx**2 - 16.d0*xx )
     . + dlnx*z2 * (  - 12.d0 - 48.d0*xx**2 + 24.d0*xx )
     . + dlnx * ( 31.d0/2.d0 + 87.d0*xx**2 - 201.d0/2.d0*xx )
     . + dlxm1**3 * (  - 35.d0/3.d0 - 70.d0/3.d0*xx**2 + 70.d0/3.d0*xx )
     . + dlxm1**2 * ( 23.d0 + 63.d0*xx**2 - 80.d0*xx )
     . + dlxm1*dli2a * ( 6.d0 + 52.d0*xx**2 - 12.d0*xx )
     . + dlxm1*z2 * ( 8.d0 + 16.d0*xx**2 - 16.d0*xx )
      partCF = partCF 
     . + dlxm1 * (  - 26.d0 - 88.d0*xx**2 + 135.d0*xx )
     . + dli2a * (  9.d0 - 40.d0*xx**2 + 24.d0*xx )
     . + dimz * (  - 8.d0 + 24.d0*xx**2 + 16.d0*xx )
     . + z2 * (  - 10.d0 + 24.d0*xx**2 - 4.d0*xx )
     . + dli3a * ( 2.d0 - 36.d0*xx**2 - 4.d0*xx )
     . + trimz * (  - 16.d0 - 32.d0*xx**2 + 32.d0*xx )
     . + s1mz * ( 22.d0 + 68.d0*xx**2 - 44.d0*xx )
     . + z3 * (  - 50.d0 - 100.d0*xx**2 + 100.d0*xx )
     . + 157.d0/4.d0 + 305.d0/4.d0*xx**2 - 221.d0/2.d0*xx

      partCF = - cf * partCF

      scaleCA = 
     .     lfac**2*( dlnx * ( 2.d0 + 8.d0*xx )
     .     + dlxm1 * ( 2.d0 + 4.d0*xx**2 - 4.d0*xx )
     .     + ( 14.d0 - 9.d0*xx**2 + 2.d0*xx + 4.d0/xx ) / 3.d0 ) +
     .     lfac*lfren*( -11.d0/3.d0 - 22.d0/3.d0*xx**2 + 22.d0/3.d0*xx )
     .     + lfac   *( dlnx**2 * (  - 4.d0 - 12.d0*xx )
     .     + dlnx*dlxm1 * ( 4.d0 - 8.d0*xx**2 + 40.d0*xx )
     .     + dlnx*dlxp1 * (  - 4.d0 - 8.d0*xx**2 - 8.d0*xx )
     .     + dlnx * ( 7.d0 + 146.d0*xx**2 + 10.d0*xx ) / 3.d0
     .     + dlxm1**2 * ( 6.d0 + 12.d0*xx**2 - 12.d0*xx )
     .     + dlxm1 * ( 40.d0/3.d0 - 98.d0/3.d0*xx**2 +
     .     64.d0/3.d0*xx + 16.d0/3.d0/xx )
     .     + dli2a * ( 12.d0 + 8.d0*xx**2 + 24.d0*xx )
     .     + dimz * (  - 4.d0 - 8.d0*xx**2 - 8.d0*xx )
     .     + z2 * (  - 8.d0 - 16.d0*xx**2 + 8.d0*xx )
     .     - 47.d0/6.d0 - 85.d0/18.d0*xx**2 + 29.d0/3.d0*xx
     .     + 44.d0/9.d0/xx  ) +
     .     lfren   *( dlnx * ( 11.d0 + 22.d0*xx**2 - 22.d0*xx ) / 3.d0
     .     + dlxm1 * ( -22.d0 - 44.d0*xx**2 + 44.d0*xx ) / 3.d0
     .     - 11.d0/6.d0 + 77.d0/6.d0*xx**2 - 11.d0*xx )

      scaleCA = CA * scaleCA

      scaleCF =
     .    lfac**2*( dlnx * (  - 3.d0 - 12.d0*xx**2 + 6.d0*xx )
     .             + dlxm1 * ( 6.d0 + 12.d0*xx**2 - 12.d0*xx )
     .             - 3.d0/2.d0 + 6.d0*xx                            ) +
     .    lfac   *( dlnx**2 * ( 4.d0 + 16.d0*xx**2 - 8.d0*xx )
     .             + dlnx*dlxm1 * (  - 20.d0 - 64.d0*xx**2 + 40.d0*xx )
     .             + dlnx * ( 5.d0 + 46.d0*xx**2 - 40.d0*xx )
     .             + dlxm1**2 * ( 18.d0 + 36.d0*xx**2 - 36.d0*xx )
     .             + dlxm1 * (  - 16.d0 - 46.d0*xx**2 + 68.d0*xx )
     .             + dli2a * (  - 24.d0*xx**2 )
     .             + z2 * (  - 4.d0 - 8.d0*xx**2 + 8.d0*xx )
     .             + 12.d0 + 11.d0*xx**2 - 34.d0*xx                  )
      
      scaleCF = CF * scaleCF

      scaleNF =
     .    lfac**2*( - 2.d0/3.d0 - 4.d0/3.d0*xx**2 + 4.d0/3.d0*xx ) +
     .    lfac*lfren*( 2.d0/3.d0 + 4.d0/3.d0*xx**2 - 4.d0/3.d0*xx ) +
     .    lfac   *( dlnx * ( 2.d0 + 4.d0*xx**2 - 4.d0*xx ) / 3.d0
     .             + dlxm1 * (  - 4.d0 - 8.d0*xx**2 + 8.d0*xx ) / 3.d0
     .             - 1.d0/3.d0 + 7.d0/3.d0*xx**2 - 2.d0*xx           ) +
     .    lfren   *( dlnx * (  - 2.d0 - 4.d0*xx**2 + 4.d0*xx ) / 3.d0
     .             + dlxm1 * ( 4.d0 + 8.d0*xx**2 - 8.d0*xx ) / 3.d0
     .             + 1.d0/3.d0 - 7.d0/3.d0*xx**2 + 2.d0*xx )

      scaleNF = NF * scaleNF

      deltabga = xx*(partCA + partCF + scaleCA + scaleCF + scaleNF)
      deltabga = deltabga/(4.d0*pi)**2

      return
      end

c.......................................................................
      double precision function deltabgf(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr
      deltabgf = 0d0
      return
      end

***
*     Regular term in gg - CGG
**
      double precision function deltagga(xx,lfh,lfr)
      implicit none
       
      double precision xx,lfh,lfr,z,lfac,lfren
      double precision xm1,xp1,dlnx,dlxm1,dlxp1
      double precision dli2a,dli3a,zcom,s1mz,dimz,trimz,smz
      double precision trimco,tripco
      double precision z2,z3,z4
      double precision delCA,delCF1,delCF2,scale
      double precision part1, part2
      double precision ca,cf,nf,pi
      double precision dilog,trilog,wgplg
      external dilog,trilog,wgplg
c     lfh here is actually -lfh,while lfr = -lfh+lfr
      lfren = -lfh+lfr
      lfac = -lfh      
      pi = 3.1415926535897932385d0
 
      ca=3.d0
      cf=4.d0/ca
      nf=5d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)
      xp1   = 1.0d0 + xx
      dlxp1 = dlog(xp1)
      dli2a = dilog(xm1)
      dli3a = trilog(xm1)
      dimz  = wgplg(1,1,-xx)
      trimz = wgplg(2,1,-xx)
      zcom  = xm1/xp1
      tripco= wgplg(2,1,zcom)
      trimco= wgplg(2,1,-zcom)
      s1mz  = wgplg(1,2,xm1)      
      smz   = wgplg(1,2,-xx)

      delCF1 = dli3a * ( 16.0d0 + 64.0d0*xx + 64.0d0*xx**2 ) +
     .       trimz * ( - 8.0d0 - 16.0d0*xx + 8.0d0*xx**2 ) +
     .       s1mz * ( - 8.0d0 - 80.0d0*xx - 56.0d0*xx**2 ) +
     .       smz * ( - 16.0d0 - 32.0d0*xx - 16.0d0*xx**2 ) +
     .       z3 * ( - 4.0d0 - 8.0d0*xx + 8.0d0*xx**2 ) +
     .       dlxm1*dli2a * ( - 16.0d0 - 64.0d0*xx - 64.0d0*xx**2 ) +
     .       dlnx*dli2a * ( - 4.0d0 - 16.0d0*xx - 16.0d0*xx**2 ) +
     .       dlxp1*dimz * ( - 16.0d0 - 32.0d0*xx - 16.0d0*xx**2 ) +
     .       dlnx*dimz * ( 16.0d0 + 32.0d0*xx + 8.0d0*xx**2 ) +
     .       dlnx*z2 * ( 12.0d0 + 40.0d0*xx + 40.0d0*xx**2 ) +
     .       dlxp1*z2 * ( - 8.0d0 - 16.0d0*xx - 8.0d0*xx**2 ) +
     .       dlxm1**2*dlnx * ( - 8.0d0 - 32.0d0*xx - 32.0d0*xx**2 ) +
     .       dlxm1*dlnx**2 * ( 4.0d0 + 16.0d0*xx + 16.0d0*xx**2 ) +
     .       dlnx**3 * ( - 2.0d0 - 16.0d0/3.0d0*xx
     .                   - 16.0d0/3.0d0*xx**2 ) +
     .       dlnx**2*dlxp1 * ( 12.0d0 + 24.0d0*xx + 12.0d0*xx**2 ) +
     .       dlnx*dlxp1**2 * ( - 8.0d0 - 16.0d0*xx - 8.0d0*xx**2 )
      delCF2 = dli2a * ( - 20.0d0 - 16.0d0*xx + 56.0d0*xx**2 ) +
     .       dimz * ( 8.0d0 + 8.0d0*xx ) +
     .       z2 * ( 20.0d0 + 36.0d0*xx - 48.0d0*xx**2 ) +
     .       dlxm1**2 * ( - 16.0d0 - 32.0d0*xx + 48.0d0*xx**2 ) +
     .       dlxm1*dlnx * ( 4.0d0 + 32.0d0*xx - 16.0d0*xx**2 ) +
     .       dlnx**2 * ( - 6.0d0 - 14.0d0*xx - 8.0d0*xx**2 ) +
     .       dlnx*dlxp1 * ( 8.0d0 + 8.0d0*xx ) +
     .       dlxm1 * ( 14.0d0 + 120.0d0*xx - 134.0d0*xx**2 ) +
     .       dlnx * ( - 23.0d0 - 64.0d0*xx + 105.0d0*xx**2 )
     .       - 32.0d0 - 66.0d0*xx + 98.0d0*xx**2

      delCA =  s1mz * ( -8.0d0*xx**2 + 16.0d0*xx - 8.0d0 ) +
     .     smz * ( 16.0d0*xx**2 + 32.0d0*xx + 16.0d0 ) +
     .     trimz * ( 24.0d0*xx**2 + 48.0d0*xx + 24.0d0 ) +
     .     z3 * ( 16.0d0*xx**2 + 32.0d0*xx + 16.0d0 ) +
     .     dimz * dlnx * ( -24.0d0*xx**2 - 48.0d0*xx - 24.0d0 ) +
     .     dimz * dlxp1 * ( 16.0d0*xx**2 + 32.0d0*xx + 16.0d0 ) +
     .     dimz * ( 16.0d0*xx**2 + 32.0d0*xx +16.0d0 )/3.0d0 +
     .     z2 * dlxp1 * ( 8.0d0*xx**2 + 16.0d0*xx + 8.0d0 ) +
     .     z2 * ( 8.0d0*xx**2 + 16.0d0*xx +8.0d0 )/3.0d0 +
     .     dlnx**2 * dlxp1 * ( -12.0d0*xx**2 - 24.0d0*xx - 12.0d0 ) +
     .     dlnx**2 * ( 50.0d0*xx**2 + 4.0d0*xx - 4.0d0 )/3.0d0 +
     .     dlnx * dlxp1**2 * ( 8.0d0*xx**2 + 16.0d0*xx + 8.0d0 ) +
     .     dlnx * dlxp1 * ( 16.0d0*xx**2 + 32.0d0*xx +
     .                       16 .0d0)/3.0d0 +
     .     dlnx * ( -50.0d0*xx**2 - 76.0d0/3.0d0*xx - 4.0d0 ) +
     .     191.0d0/3.0d0*xx**2 - 48.0d0*xx - 47.0d0/3.0d0 
      delCA =  ca**2 / ( ca**2 - 1.0d0 ) * delCA

      scale = ( dlnx * ( - 2.0d0 - 8.0d0*xx - 8.0d0*xx**2 )
     .           - 4.0d0 - 8.0d0*xx + 12.0d0*xx**2  ) * lfac**2 +
     .         ( dli2a * ( - 8.0d0 - 32.0d0*xx - 32.0d0*xx**2 ) +
     .           dlxm1*dlnx * ( - 8.0d0 - 32.0d0*xx - 32.0d0*xx**2 ) +
     .           dlnx**2 * ( 2.0d0 + 8.0d0*xx + 8.0d0*xx**2 ) +
     .           dlxm1 * ( - 16.0d0 - 32.0d0*xx + 48.0d0*xx**2 ) +
     .           dlnx * ( 2.0d0 + 16.0d0*xx - 8.0d0*xx**2 ) +
     .           7.0d0 + 60.0d0*xx - 67.0d0*xx**2) * lfac

      deltagga = xx*(delCF1 + delCF2 + delCA + scale)
      deltagga = deltagga/(4.d0*pi)**2

      return
      end

c.......................................................................
      double precision function deltaggf(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr
      deltaggf = 0.d0
      return
      end

***
*     Regular term in qq - CQQ
**
      double precision function CQQIV(xx,lfh,lfr)
      implicit none
      double precision partIV
      double precision xx,lfh,lfr,z
      double precision xm1,xp1,dlnx,dlxm1,dlxp1
      double precision dli2a,dli3a,zcom,s1mz,dimz,trimz,smz
      double precision trimco,tripco
      double precision part1,part2,scale
      double precision z2,z3,z4
      double precision ca,cf,nf
      double precision dilog,trilog,wgplg
      double precision cax,cv,cd,sw
      external dilog,trilog,wgplg

      sw  = 0.2228972225239183d0
      cax = 1.d0
      cv  = -1.d0 + 4.d0/3.d0*sw
      cd  = 1.d0 + cv**2
      
      ca=3.d0
      cf=4.d0/ca
      nf=5d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)
      xp1   = 1.0d0 + xx
      dlxp1 = dlog(xp1)
      dli2a = dilog(xm1)
      dli3a = trilog(xm1)
      dimz  = wgplg(1,1,-xx)
      trimz = wgplg(2,1,-xx)
      zcom  = xm1/xp1
      tripco= wgplg(2,1,zcom)
      trimco= wgplg(2,1,-zcom)
      s1mz  = wgplg(1,2,xm1)
      smz   = wgplg(1,2,-xx)

      Z = xx

      
      partIV =
     .   dlnx**3 * ( 4.D0/3.D0*Z )
     . + dlnx**2*dlxp1 * (  - 20.D0 - 10.D0*Z - 20.D0/Z )
     . + dlnx**2 * ( 13.D0*Z )
     . + dlnx*dlxp1**2 * ( 24.D0 + 12.D0*Z + 24.D0/Z )
     . + dlnx*dlxp1 * (  - 20.D0 - 20.D0*Z )
     . + dlnx*dli2a * (  - 20.D0 + 2.D0*Z )
     . + dlnx*DIMZ * (  - 16.D0*Z - 40.D0/Z )
     . + dlnx*z2 * (  - 20.D0 - 2.D0*Z )
     . + dlnx * ( 20.D0 + 16.D0*Z )
     . + dlxp1*DIMZ * ( 48.D0 + 24.D0*Z + 48.D0/Z )
     . + dlxp1*z2 * ( 24.D0 + 12.D0*Z + 24.D0/Z )
      partIV = partIV
     . + dli2a * (  - 10.D0 + 8.D0*Z )
     . + DIMZ * (  - 20.D0 - 20.D0*Z )
     . + z2 * (  - 10.D0 - 10.D0*Z )
     . + dli3a * ( 12.D0 - 6.D0*Z - 8.D0/Z )
     . + TRIMZ * (  - 40.D0 + 12.D0*Z + 40.D0/Z )
     . + S1MZ * (  - 16.D0 - 8.D0*Z - 16.D0/Z )
     . + SMZ * ( 48.D0 + 24.D0*Z + 48.D0/Z )
     . + z3 * (  - 36.D0 + 6.D0*Z + 24.D0/Z )
     . + 40.D0 - 40.D0*Z
      partIV = cf * partIV
      CQQIV  = partIV
      return
      end


      double precision function CQQIA(xx,lfh,lfr)
      implicit none
      double precision partIA
      double precision xx,lfh,lfr,z
      double precision xm1,xp1,dlnx,dlxm1,dlxp1
      double precision dli2a,dli3a,zcom,s1mz,dimz,trimz,smz
      double precision trimco,tripco
      double precision part1,part2,scale
      double precision z2,z3,z4
      double precision ca,cf,nf
      double precision dilog,trilog,wgplg
      double precision cax,cv,cd,sw
      external dilog,trilog,wgplg

      sw  = 0.2228972225239183d0
      cax = 1.d0
      cv  = -1.d0 + 4.d0/3.d0*sw
      cd  = 1.d0 + cv**2

      
      ca=3.d0
      cf=4.d0/ca
      nf=5d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)
      xp1   = 1.0d0 + xx
      dlxp1 = dlog(xp1)
      dli2a = dilog(xm1)
      dli3a = trilog(xm1)
      dimz  = wgplg(1,1,-xx)
      trimz = wgplg(2,1,-xx)
      zcom  = xm1/xp1
      tripco= wgplg(2,1,zcom)
      trimco= wgplg(2,1,-zcom)
      s1mz  = wgplg(1,2,xm1)
      smz   = wgplg(1,2,-xx)

      Z = xx 

c---  CQQIA

      partIA =
     .   dlnx**3 * (  - 4.D0/3.D0*Z )
     . + dlnx**2*dlxp1 * ( 20.D0 + 10.D0*Z )
     . + dlnx**2 * (  - Z )
     . + dlnx*dlxp1**2 * (  - 24.D0 - 12.D0*Z )
     . + dlnx*dlxp1 * ( 4.D0 + 4.D0*Z )
     . + dlnx*dli2a * ( 4.D0 + 6.D0*Z )
     . + dlnx*DIMZ * ( 32.D0 )
     . + dlnx*z2 * ( 4.D0 + 10.D0*Z )  - dlnx * 4.D0
     . + dlxp1*DIMZ * (  - 48.D0 - 24.D0*Z )
     . + dlxp1*z2 * (  - 24.D0 - 12.D0*Z )
     . + dli2a * ( 2.D0 )
     . + DIMZ * ( 4.D0 + 4.D0*Z )
     . + z2 * ( 2.D0 + 2.D0*Z )
     . + dli3a * ( 4.D0 - 2.D0*Z )
     . + TRIMZ * (  - 24.D0 + 20.D0*Z )
     . + S1MZ * ( 16.D0 + 8.D0*Z )
     . + SMZ * (  - 48.D0 - 24.D0*Z )
     . + z3 * (  - 12.D0 + 18.D0*Z ) - 8.D0 + 8.D0*Z
      partIA = cf * partIA
      CQQIA = partIA
      return
      end

      

      double precision function CQQ(xx,lfh,lfr)
      implicit none
      
      double precision xx,lfh,lfr,z
      double precision xm1,xp1,dlnx,dlxm1,dlxp1
      double precision dli2a,dli3a,zcom,s1mz,dimz,trimz,smz
      double precision trimco,tripco
      double precision part1,part2,scale
      double precision z2,z3,z4
      double precision ca,cf,nf
      double precision dilog,trilog,wgplg
      double precision cax,cv,cd,sw
      external dilog,trilog,wgplg

      sw  = 0.2228972225239183d0
      cax = 1.d0
      cv  = -1.d0 + 4.d0/3.d0*sw
      cd  = 1.d0 + cv**2

      
      ca=3.d0
      cf=4.d0/ca
      nf=5d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)
      xp1   = 1.0d0 + xx
      dlxp1 = dlog(xp1)
      dli2a = dilog(xm1)
      dli3a = trilog(xm1)
      dimz  = wgplg(1,1,-xx)
      trimz = wgplg(2,1,-xx)
      zcom  = xm1/xp1
      tripco= wgplg(2,1,zcom)
      trimco= wgplg(2,1,-zcom)
      s1mz  = wgplg(1,2,xm1)
      smz   = wgplg(1,2,-xx)

      Z = xx

      part1=-8.0d0*dli3a + 12.0d0*s1mz + 8.0d0*dlxm1*dli2a +
     .       2.0d0*dlnx*dli2a - 4.0d0*dlnx*z2 + 1.5d0*dlnx**3
     .      -4.0d0*dlnx**2*dlxm1 + 4.0d0*dlnx*dlxm1**2
      part1= xp1 * part1
      part2=dli2a * ( 13.0d0 + 8.0d0/3.0d0*z**2 + 5.0d0*z +
     .                16.0d0/3.0d0/z )
     .    + z2 * ( - 2.0d0 + 8.0d0/3.0d0*z**2 + 2.0d0*z
     .                - 8.0d0/3.0d0/z )
     .    + dlxm1**2 * ( 2.0d0 - 8.0d0/3.0d0*z**2 - 2.0d0*z +
     .                    8.0d0/3.0d0/z )
     .    + dlnx*dlxm1 * ( 6.0d0 + 8.0d0*z**2 + 12.0d0*z )
     .    + dlnx**2 * ( - 5.0d0/4.0d0 - 10.0d0/3.0d0*z**2
     .                  - 25.0d0/4.0d0*z )
     .    + dlxm1 * ( - 26.0d0/3.0d0 - 44.0d0/9.0d0*z**2 +
     .                   26.0d0/3.0d0*z + 44.0d0/9.0d0/z )
     .    + dlnx * ( 115.0d0/6.0d0 + 10.0d0/9.0d0*z**2
     .             - 8.0d0/3.0d0*z )
     .    + 593.0d0/36.0d0 + 703.0d0/108.0d0*z**2
     .    - 433.0d0/18.0d0*z + 29.0d0/27.0d0/z
      

      scale =lfh**2 * ( dlnx * ( 1.0d0 + z ) +
     .     0.5d0 - 2.0d0/3.0d0*z**2 - 0.5d0*z +
     .     2.0d0/3.0d0/z ) +
     .     lfh * ( dli2a * ( 4.0d0 + 4.0d0*z ) +
     .     dlnx**2 * ( - 2.0d0 - 2.0d0*z ) +
     .     dlnx*dlxm1 * ( 4.0d0 + 4.0d0*z ) +
     .     dlnx * ( 3.0d0 + 4.0d0*z**2 + 6.0d0*z ) +
     .     dlxm1 * (  2.0d0 - 8.0d0/3.0d0*z**2
     .     - 2.0d0*z + 8.0d0/3.0d0/z )
     .     -13.0d0/3.0d0 - 22.0d0/9.0d0*z**2 +
     .     13.0d0/3.0d0*z + 22.0d0/9.0d0/z )
      
      CQQ = ( ca**2 - 1.0d0 )/ ca *( part1 + part2 + scale )
      
      return
      end

      DOUBLE PRECISION FUNCTION CQBBB(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr,z
      double precision logz,s1mz,dimz,trimz,smz
      double precision log1pz,di1mz
      double precision trimco,tripco
      double precision part1,part2,scale
      double precision z2,z3,z4
      double precision ca,cf,nf
      double precision dilog,trilog,wgplg
      external dilog,trilog,wgplg
      
      ca=3.d0
      cf=4.d0/ca
      nf=5d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      Z = xx

      LOGZ=DLOG(Z)
      LOG1PZ=DLOG(1.0D0+Z)
      DIMZ=WGPLG(1,1,-Z)
      CQBBB =CF * (
     .   LOGZ**2 * ( 4.D0/3.D0 + 4.D0/3.D0*Z**2 + 8.D0/3.D0*Z )
     . + LOGZ*LOG1PZ * (  - 16.D0/3.D0 - 16.D0/3.D0*Z**2
     .                    - 32.D0/3.D0*Z )
     . + LOGZ * ( 4.D0 + 4.D0*Z**2 + 16.D0/3.D0*Z )
     . + DIMZ * (  - 16.D0/3.D0 - 16.D0/3.D0*Z**2 - 32.D0/3.D0*Z )
     . + Z2 * (  - 8.D0/3.D0 - 8.D0/3.D0*Z**2 - 16.D0/3.D0*Z )
     . + 20.D0/3.D0 - 20.D0/3.D0*Z**2 )
      RETURN
      END

      DOUBLE PRECISION FUNCTION CQBBCD(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr,z
      double precision logz,s1mz,dimz,trimz,smz
      double precision log1pz,di1mz
      double precision trimco,tripco
      double precision part1,part2,scale
      double precision z2,z3,z4
      double precision ca,cf,nf
      double precision dilog,trilog,wgplg
      external dilog,trilog,wgplg
      
      ca=3.d0
      cf=4.d0/ca
      nf=5d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      Z = xx

      LOGZ=DLOG(Z)
      LOG1PZ=DLOG(1.0D0+Z)
      DIMZ=WGPLG(1,1,-Z)
      DI1MZ=WGPLG(1,1,1.0D0-Z)
      TRIMZ=WGPLG(2,1,-Z)
      SMZ=WGPLG(1,2,-Z)
      S1MZ=WGPLG(1,2,1.0D0-Z)
      CQBBCD =CF * ( CF - CA/2 ) * (
     .   LOGZ**3 * ( 4.D0/3.D0 + 4.D0/3.D0*Z**2 + 16.D0/3.D0*Z )
     . + LOGZ**2*LOG1PZ * ( 20.D0 + 20.D0*Z**2 + 40.D0*Z )
     . + LOGZ**2 * ( 12.D0 - 30.D0*Z**2 - 16.D0*Z )
     . + LOGZ*LOG1PZ**2 * (  - 24.D0 - 24.D0*Z**2 - 48.D0*Z )
     . + LOGZ*LOG1PZ * ( 24.D0 + 24.D0*Z**2 + 48.D0*Z )
     . + LOGZ*DI1MZ * ( 16.D0 + 16.D0*Z**2 + 48.D0*Z )
     . + LOGZ*DIMZ * ( 24.D0 + 24.D0*Z**2 + 48.D0*Z )
     . + LOGZ*Z2 * ( 8.D0 + 8.D0*Z**2 + 16.D0*Z )
     . + LOGZ * ( 36.D0 + 44.D0*Z )
     . + LOG1PZ*DIMZ * (  - 48.D0 - 48.D0*Z**2 - 96.D0*Z )
     . + LOG1PZ*Z2 * (  - 24.D0 - 24.D0*Z**2 - 48.D0*Z )
     . + DI1MZ * ( 36.D0 - 36.D0*Z**2 )
     . + DIMZ * ( 24.D0 + 24.D0*Z**2 + 48.D0*Z )
     . + Z2 * ( 12.D0 + 12.D0*Z**2 + 24.D0*Z )
     . + TRIMZ * (  - 8.D0 - 8.D0*Z**2 - 16.D0*Z )
     . + S1MZ * ( 32.D0 + 32.D0*Z**2 + 96.D0*Z )
     . + SMZ * (  - 48.D0 - 48.D0*Z**2 - 96.D0*Z )
     . + 54.D0 - 26.D0*Z**2 - 28.D0*Z )
      RETURN
      END

      double precision function CQBAB(xx,lfh,lfr)
      implicit none
      double precision partBA
      double precision xx,lfh,lfr,z
      double precision xm1,xp1,dlnx,dlxm1,dlxp1
      double precision dli2a,dli3a,zcom,s1mz,dimz,trimz,smz
      double precision trimco,tripco
      double precision part1,part2,scale
      double precision z2,z3,z4
      double precision ca,cf,nf
      double precision dilog,trilog,wgplg
      double precision cax,cv,cd,sw
      external dilog,trilog,wgplg
      ca=3.d0
      cf=4.d0/ca
      nf=5d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)
      xp1   = 1.0d0 + xx
      dlxp1 = dlog(xp1)
      dli2a = dilog(xm1)
      dli3a = trilog(xm1)
      dimz  = wgplg(1,1,-xx)
      trimz = wgplg(2,1,-xx)
      zcom  = xm1/xp1
      tripco= wgplg(2,1,zcom)
      trimco= wgplg(2,1,-zcom)
      s1mz  = wgplg(1,2,xm1)
      smz   = wgplg(1,2,-xx)

      Z = xx

      
      partBA =  dlnx * ( - 8.d0 + 8.d0*xx + 16.d0/xm1 )
     . + 24.d0 - 8.d0*xx  

      partBA =CF * partBA

      CQBAB = partBA
      
      return
      end
      
      double precision function deltabba(xx,lfh,lfr)
      implicit none
      
      double precision xx,lfh,lfr,z,lfac,lfren,pi
      double precision CQQIV,CQQ,CQQIA,CQBBB
      double precision CQBBCD,CQQCRF,CQQCRI
      double precision cax,cv,cd,sw
      external CQQIV,CQQ,CQQIA,CQBBB
      external CQBBCD,CQQCRF,CQQCRI

      sw  = 0.2228972225239183d0
      cax = 1.d0
      cv  = -1.d0 + 4.d0/3.d0*sw
      cd  = 1.d0 + cv**2
c     lfh here is actually -lfh,while lfr = -lfh+lfr
      lfren = -lfh+lfr
      lfac = -lfh      
      pi = 3.1415926535897932385d0
 
      deltabba = xx*(cd*CQBBB(xx,lfac,lfren)
     ?     + 2.d0*(cd*(CQQ(xx,lfac,lfren)
     ?     + CQQCRF(xx,lfac,lfren)
     ?     + CQQCRI(xx,lfac,lfren))
     ?     +cax**2*CQQIA(xx,lfac,lfren)
     ?     +cv**2 *CQQIV(xx,lfac,lfren)))
      deltabba = deltabba/cd
      deltabba = deltabba/(4.d0*pi)**2
      return
      end

c.......................................................................
      double precision function deltabbf(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr,CQQ
      external CQQ
      deltabbf = 0.d0
      return
      end

*******************************************

      double precision function deltabqa(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr,lfac,lfren,pi
      double precision CQQ
      external CQQ
c     lfh here is actually -lfh,while lfr = -lfh+lfr
      lfren = -lfh+lfr
      lfac = -lfh      
      pi = 3.1415926535897932385d0
      deltabqa = xx*(CQQ(xx,lfac,lfren) )
      deltabqa = deltabqa/(4.d0*pi)**2
      end
c.......................................................................

      double precision function deltabqf(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr

      deltabqf = 0.d0
      end
      
*******************************************

      double precision function deltaqqba(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr,lfac,lfren,pi
      double precision CQBBB
      external CQBBB
c     lfh here is actually -lfh,while lfr = -lfh+lfr
      lfren = -lfh+lfr
      lfac = -lfh      
      pi = 3.1415926535897932385d0
 
      deltaqqba = xx*CQBBB(xx,lfac,lfren)
      deltaqqba = deltaqqba/(4.0d0*pi)**2
      return
      end

c.......................................................................
      double precision function deltaqqbf(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr
      deltaqqbf = 0.d0
      return
      end
c......................................................................
c......................................................................
c................... FUNCTIONS ........................................
c......................................................................
c......................................................................
      double precision function dilog(xx)
c..
c..   Dilogarithm: dilog(x) = Li_2(x)
c..
      double precision xx
      double precision wgplg
      dilog = wgplg(1,1,xx)
      end

c..   ------------------------------------------------------------

C-}}}
C-{{{ function trilog:

c..   ------------------------------------------------------------

      double precision function trilog(xx)
c..
c..   Trilogarithm: trilog(x) = Li_3(x)
c..
      double precision xx
      double precision wgplg
      trilog = wgplg(2,1,xx)
      end

c..   ------------------------------------------------------------

C-}}}
C-{{{ FUNCTION WGPLG:

      double precision FUNCTION WGPLG(N,P,X)
c..
c..   Nielson's Generalized Polylogarithm
c..
c..   taken from zwprod.f (R. Hamberg, T. Matsuura and W.L. van Neerven)
c..   (originally from CERNLIB?)
c..
      INTEGER P,P1,NC(10),INDEX(31)
      DOUBLE PRECISION FCT(0:4),SGN(0:4),U(0:4),S1(4,4),C(4,4)
      DOUBLE PRECISION A(0:30,10)
      DOUBLE PRECISION X,X1,H,ALFA,R,Q,C1,C2,B0,B1,B2,ZERO,HALF

      DOUBLE PRECISION V(0:5),SK,SM

      DATA FCT /1.0D0,1.0D0,2.0D0,6.0D0,24.0D0/
      DATA SGN /1.0D0,-1.0D0,1.0D0,-1.0D0,1.0D0/
      DATA ZERO /0.0D0/, HALF /0.5D0/
      DATA C1 /1.33333 33333 333D0/, C2 /0.33333 33333 3333D0/

      DATA S1(1,1) /1.64493 40668 482D0/
      DATA S1(1,2) /1.20205 69031 596D0/
      DATA S1(1,3) /1.08232 32337 111D0/
      DATA S1(1,4) /1.03692 77551 434D0/
      DATA S1(2,1) /1.20205 69031 596D0/
      DATA S1(2,2) /2.70580 80842 778D-1/
      DATA S1(2,3) /9.65511 59989 444D-2/
      DATA S1(3,1) /1.08232 32337 111D0/
      DATA S1(3,2) /9.65511 59989 444D-2/
      DATA S1(4,1) /1.03692 77551 434D0/

      DATA C(1,1) / 1.64493 40668 482D0/
      DATA C(1,2) / 1.20205 69031 596D0/
      DATA C(1,3) / 1.08232 32337 111D0/
      DATA C(1,4) / 1.03692 77551 434D0/
      DATA C(2,1) / 0.00000 00000 000D0/
      DATA C(2,2) /-1.89406 56589 945D0/
      DATA C(2,3) /-3.01423 21054 407D0/
      DATA C(3,1) / 1.89406 56589 945D0/
      DATA C(3,2) / 3.01423 21054 407D0/
      DATA C(4,1) / 0.00000 00000 000D0/

      DATA INDEX /1,2,3,4,6*0,5,6,7,7*0,8,9,8*0,10/

      DATA NC /24,26,28,30,22,24,26,19,22,17/

      DATA A( 0,1) / .96753 21504 3498D0/
      DATA A( 1,1) / .16607 30329 2785D0/
      DATA A( 2,1) / .02487 93229 2423D0/
      DATA A( 3,1) / .00468 63619 5945D0/
      DATA A( 4,1) / .00100 16274 9616D0/
      DATA A( 5,1) / .00023 20021 9609D0/
      DATA A( 6,1) / .00005 68178 2272D0/
      DATA A( 7,1) / .00001 44963 0056D0/
      DATA A( 8,1) / .00000 38163 2946D0/
      DATA A( 9,1) / .00000 10299 0426D0/
      DATA A(10,1) / .00000 02835 7538D0/
      DATA A(11,1) / .00000 00793 8705D0/
      DATA A(12,1) / .00000 00225 3670D0/
      DATA A(13,1) / .00000 00064 7434D0/
      DATA A(14,1) / .00000 00018 7912D0/
      DATA A(15,1) / .00000 00005 5029D0/
      DATA A(16,1) / .00000 00001 6242D0/
      DATA A(17,1) / .00000 00000 4827D0/
      DATA A(18,1) / .00000 00000 1444D0/
      DATA A(19,1) / .00000 00000 0434D0/
      DATA A(20,1) / .00000 00000 0131D0/
      DATA A(21,1) / .00000 00000 0040D0/
      DATA A(22,1) / .00000 00000 0012D0/
      DATA A(23,1) / .00000 00000 0004D0/
      DATA A(24,1) / .00000 00000 0001D0/

      DATA A( 0,2) / .95180 88912 7832D0/
      DATA A( 1,2) / .43131 13184 6532D0/
      DATA A( 2,2) / .10002 25071 4905D0/
      DATA A( 3,2) / .02442 41559 5220D0/
      DATA A( 4,2) / .00622 51246 3724D0/
      DATA A( 5,2) / .00164 07883 1235D0/
      DATA A( 6,2) / .00044 40792 0265D0/
      DATA A( 7,2) / .00012 27749 4168D0/
      DATA A( 8,2) / .00003 45398 1284D0/
      DATA A( 9,2) / .00000 98586 9565D0/
      DATA A(10,2) / .00000 28485 6995D0/
      DATA A(11,2) / .00000 08317 0847D0/
      DATA A(12,2) / .00000 02450 3950D0/
      DATA A(13,2) / .00000 00727 6496D0/
      DATA A(14,2) / .00000 00217 5802D0/
      DATA A(15,2) / .00000 00065 4616D0/
      DATA A(16,2) / .00000 00019 8033D0/
      DATA A(17,2) / .00000 00006 0204D0/
      DATA A(18,2) / .00000 00001 8385D0/
      DATA A(19,2) / .00000 00000 5637D0/
      DATA A(20,2) / .00000 00000 1735D0/
      DATA A(21,2) / .00000 00000 0536D0/
      DATA A(22,2) / .00000 00000 0166D0/
      DATA A(23,2) / .00000 00000 0052D0/
      DATA A(24,2) / .00000 00000 0016D0/
      DATA A(25,2) / .00000 00000 0005D0/
      DATA A(26,2) / .00000 00000 0002D0/

      DATA A( 0,3) / .98161 02799 1365D0/
      DATA A( 1,3) / .72926 80632 0726D0/
      DATA A( 2,3) / .22774 71490 9321D0/
      DATA A( 3,3) / .06809 08329 6197D0/
      DATA A( 4,3) / .02013 70118 3064D0/
      DATA A( 5,3) / .00595 47848 0197D0/
      DATA A( 6,3) / .00176 76901 3959D0/
      DATA A( 7,3) / .00052 74821 8502D0/
      DATA A( 8,3) / .00015 82746 1460D0/
      DATA A( 9,3) / .00004 77492 2076D0/
      DATA A(10,3) / .00001 44792 0408D0/
      DATA A(11,3) / .00000 44115 4886D0/
      DATA A(12,3) / .00000 13500 3870D0/
      DATA A(13,3) / .00000 04148 1779D0/
      DATA A(14,3) / .00000 01279 3307D0/
      DATA A(15,3) / .00000 00395 9070D0/
      DATA A(16,3) / .00000 00122 9055D0/
      DATA A(17,3) / .00000 00038 2658D0/
      DATA A(18,3) / .00000 00011 9459D0/
      DATA A(19,3) / .00000 00003 7386D0/
      DATA A(20,3) / .00000 00001 1727D0/
      DATA A(21,3) / .00000 00000 3687D0/
      DATA A(22,3) / .00000 00000 1161D0/
      DATA A(23,3) / .00000 00000 0366D0/
      DATA A(24,3) / .00000 00000 0116D0/
      DATA A(25,3) / .00000 00000 0037D0/
      DATA A(26,3) / .00000 00000 0012D0/
      DATA A(27,3) / .00000 00000 0004D0/
      DATA A(28,3) / .00000 00000 0001D0/

      DATA A( 0,4) /1.06405 21184 614 D0/
      DATA A( 1,4) /1.06917 20744 981 D0/
      DATA A( 2,4) / .41527 19325 1768D0/
      DATA A( 3,4) / .14610 33293 6222D0/
      DATA A( 4,4) / .04904 73264 8784D0/
      DATA A( 5,4) / .01606 34086 0396D0/
      DATA A( 6,4) / .00518 88935 0790D0/
      DATA A( 7,4) / .00166 29871 7324D0/
      DATA A( 8,4) / .00053 05827 9969D0/
      DATA A( 9,4) / .00016 88702 9251D0/
      DATA A(10,4) / .00005 36832 8059D0/
      DATA A(11,4) / .00001 70592 3313D0/
      DATA A(12,4) / .00000 54217 4374D0/
      DATA A(13,4) / .00000 17239 4082D0/
      DATA A(14,4) / .00000 05485 3275D0/
      DATA A(15,4) / .00000 01746 7795D0/
      DATA A(16,4) / .00000 00556 7550D0/
      DATA A(17,4) / .00000 00177 6234D0/
      DATA A(18,4) / .00000 00056 7224D0/
      DATA A(19,4) / .00000 00018 1313D0/
      DATA A(20,4) / .00000 00005 8012D0/
      DATA A(21,4) / .00000 00001 8579D0/
      DATA A(22,4) / .00000 00000 5955D0/
      DATA A(23,4) / .00000 00000 1911D0/
      DATA A(24,4) / .00000 00000 0614D0/
      DATA A(25,4) / .00000 00000 0197D0/
      DATA A(26,4) / .00000 00000 0063D0/
      DATA A(27,4) / .00000 00000 0020D0/
      DATA A(28,4) / .00000 00000 0007D0/
      DATA A(29,4) / .00000 00000 0002D0/
      DATA A(30,4) / .00000 00000 0001D0/

      DATA A( 0,5) / .97920 86066 9175D0/
      DATA A( 1,5) / .08518 81314 8683D0/
      DATA A( 2,5) / .00855 98522 2013D0/
      DATA A( 3,5) / .00121 17721 4413D0/
      DATA A( 4,5) / .00020 72276 8531D0/
      DATA A( 5,5) / .00003 99695 8691D0/
      DATA A( 6,5) / .00000 83806 4065D0/
      DATA A( 7,5) / .00000 18684 8945D0/
      DATA A( 8,5) / .00000 04366 6087D0/
      DATA A( 9,5) / .00000 01059 1733D0/
      DATA A(10,5) / .00000 00264 7892D0/
      DATA A(11,5) / .00000 00067 8700D0/
      DATA A(12,5) / .00000 00017 7654D0/
      DATA A(13,5) / .00000 00004 7342D0/
      DATA A(14,5) / .00000 00001 2812D0/
      DATA A(15,5) / .00000 00000 3514D0/
      DATA A(16,5) / .00000 00000 0975D0/
      DATA A(17,5) / .00000 00000 0274D0/
      DATA A(18,5) / .00000 00000 0077D0/
      DATA A(19,5) / .00000 00000 0022D0/
      DATA A(20,5) / .00000 00000 0006D0/
      DATA A(21,5) / .00000 00000 0002D0/
      DATA A(22,5) / .00000 00000 0001D0/

      DATA A( 0,6) / .95021 85196 3952D0/
      DATA A( 1,6) / .29052 52916 1433D0/
      DATA A( 2,6) / .05081 77406 1716D0/
      DATA A( 3,6) / .00995 54376 7280D0/
      DATA A( 4,6) / .00211 73389 5031D0/
      DATA A( 5,6) / .00047 85947 0550D0/
      DATA A( 6,6) / .00011 33432 1308D0/
      DATA A( 7,6) / .00002 78473 3104D0/
      DATA A( 8,6) / .00000 70478 8108D0/
      DATA A( 9,6) / .00000 18278 8740D0/
      DATA A(10,6) / .00000 04838 7492D0/
      DATA A(11,6) / .00000 01303 3842D0/
      DATA A(12,6) / .00000 00356 3769D0/
      DATA A(13,6) / .00000 00098 7174D0/
      DATA A(14,6) / .00000 00027 6586D0/
      DATA A(15,6) / .00000 00007 8279D0/
      DATA A(16,6) / .00000 00002 2354D0/
      DATA A(17,6) / .00000 00000 6435D0/
      DATA A(18,6) / .00000 00000 1866D0/
      DATA A(19,6) / .00000 00000 0545D0/
      DATA A(20,6) / .00000 00000 0160D0/
      DATA A(21,6) / .00000 00000 0047D0/
      DATA A(22,6) / .00000 00000 0014D0/
      DATA A(23,6) / .00000 00000 0004D0/
      DATA A(24,6) / .00000 00000 0001D0/

      DATA A( 0,7) / .95064 03218 6777D0/
      DATA A( 1,7) / .54138 28546 5171D0/
      DATA A( 2,7) / .13649 97959 0321D0/
      DATA A( 3,7) / .03417 94232 8207D0/
      DATA A( 4,7) / .00869 02788 3583D0/
      DATA A( 5,7) / .00225 28408 4155D0/
      DATA A( 6,7) / .00059 51608 9806D0/
      DATA A( 7,7) / .00015 99561 7766D0/
      DATA A( 8,7) / .00004 36521 3096D0/
      DATA A( 9,7) / .00001 20747 4688D0/
      DATA A(10,7) / .00000 33801 8176D0/
      DATA A(11,7) / .00000 09563 2476D0/
      DATA A(12,7) / .00000 02731 3129D0/
      DATA A(13,7) / .00000 00786 6968D0/
      DATA A(14,7) / .00000 00228 3195D0/
      DATA A(15,7) / .00000 00066 7205D0/
      DATA A(16,7) / .00000 00019 6191D0/
      DATA A(17,7) / .00000 00005 8018D0/
      DATA A(18,7) / .00000 00001 7246D0/
      DATA A(19,7) / .00000 00000 5151D0/
      DATA A(20,7) / .00000 00000 1545D0/
      DATA A(21,7) / .00000 00000 0465D0/
      DATA A(22,7) / .00000 00000 0141D0/
      DATA A(23,7) / .00000 00000 0043D0/
      DATA A(24,7) / .00000 00000 0013D0/
      DATA A(25,7) / .00000 00000 0004D0/
      DATA A(26,7) / .00000 00000 0001D0/

      DATA A( 0,8) / .98800 01167 2229D0/
      DATA A( 1,8) / .04364 06760 9601D0/
      DATA A( 2,8) / .00295 09117 8278D0/
      DATA A( 3,8) / .00031 47780 9720D0/
      DATA A( 4,8) / .00004 31484 6029D0/
      DATA A( 5,8) / .00000 69381 8230D0/
      DATA A( 6,8) / .00000 12464 0350D0/
      DATA A( 7,8) / .00000 02429 3628D0/
      DATA A( 8,8) / .00000 00504 0827D0/
      DATA A( 9,8) / .00000 00109 9075D0/
      DATA A(10,8) / .00000 00024 9467D0/
      DATA A(11,8) / .00000 00005 8540D0/
      DATA A(12,8) / .00000 00001 4127D0/
      DATA A(13,8) / .00000 00000 3492D0/
      DATA A(14,8) / .00000 00000 0881D0/
      DATA A(15,8) / .00000 00000 0226D0/
      DATA A(16,8) / .00000 00000 0059D0/
      DATA A(17,8) / .00000 00000 0016D0/
      DATA A(18,8) / .00000 00000 0004D0/
      DATA A(19,8) / .00000 00000 0001D0/

      DATA A( 0,9) / .95768 50654 6350D0/
      DATA A( 1,9) / .19725 24967 9534D0/
      DATA A( 2,9) / .02603 37031 3918D0/
      DATA A( 3,9) / .00409 38216 8261D0/
      DATA A( 4,9) / .00072 68170 7110D0/
      DATA A( 5,9) / .00014 09187 9261D0/
      DATA A( 6,9) / .00002 92045 8914D0/
      DATA A( 7,9) / .00000 63763 1144D0/
      DATA A( 8,9) / .00000 14516 7850D0/
      DATA A( 9,9) / .00000 03420 5281D0/
      DATA A(10,9) / .00000 00829 4302D0/
      DATA A(11,9) / .00000 00206 0784D0/
      DATA A(12,9) / .00000 00052 2823D0/
      DATA A(13,9) / .00000 00013 5066D0/
      DATA A(14,9) / .00000 00003 5451D0/
      DATA A(15,9) / .00000 00000 9436D0/
      DATA A(16,9) / .00000 00000 2543D0/
      DATA A(17,9) / .00000 00000 0693D0/
      DATA A(18,9) / .00000 00000 0191D0/
      DATA A(19,9) / .00000 00000 0053D0/
      DATA A(20,9) / .00000 00000 0015D0/
      DATA A(21,9) / .00000 00000 0004D0/
      DATA A(22,9) / .00000 00000 0001D0/

      DATA A( 0,10) / .99343 65167 1347D0/
      DATA A( 1,10) / .02225 77012 6826D0/
      DATA A( 2,10) / .00101 47557 4703D0/
      DATA A( 3,10) / .00008 17515 6250D0/
      DATA A( 4,10) / .00000 89997 3547D0/
      DATA A( 5,10) / .00000 12082 3987D0/
      DATA A( 6,10) / .00000 01861 6913D0/
      DATA A( 7,10) / .00000 00317 4723D0/
      DATA A( 8,10) / .00000 00058 5215D0/
      DATA A( 9,10) / .00000 00011 4739D0/
      DATA A(10,10) / .00000 00002 3652D0/
      DATA A(11,10) / .00000 00000 5082D0/
      DATA A(12,10) / .00000 00000 1131D0/
      DATA A(13,10) / .00000 00000 0259D0/
      DATA A(14,10) / .00000 00000 0061D0/
      DATA A(15,10) / .00000 00000 0015D0/
      DATA A(16,10) / .00000 00000 0004D0/
      DATA A(17,10) / .00000 00000 0001D0/

      IF(N .LT. 1 .OR. N .GT. 4 .OR. P .LT. 1 .OR. P .GT. 4 .OR.
     1   N+P .GT. 5) THEN
       WGPLG=ZERO
       PRINT 1000, N,P
       RETURN
      END IF
      IF(X .EQ. SGN(0)) THEN
       WGPLG=S1(N,P)
       RETURN
      END IF

      IF(X .GT. FCT(2) .OR. X .LT. SGN(1)) THEN
       X1=SGN(0)/X
       H=C1*X1+C2
       ALFA=H+H
       V(0)=SGN(0)
       V(1)=LOG(DCMPLX(-X,ZERO))
       DO 33 L = 2,N+P
   33  V(L)=V(1)*V(L-1)/L
       SK=ZERO
       DO 34 K = 0,P-1
       P1=P-K
       R=X1**P1/(FCT(P1)*FCT(N-1))
       SM=ZERO
       DO 35 M = 0,K
       N1=N+K-M
       L=INDEX(10*N1+P1-10)
       B1=ZERO
       B2=ZERO
       DO 31 I = NC(L),0,-1
       B0=A(I,L)+ALFA*B1-B2
       B2=B1
   31  B1=B0
       Q=(FCT(N1-1)/FCT(K-M))*(B0-H*B2)*R/P1**N1
   35  SM=SM+V(M)*Q
   34  SK=SK+SGN(K)*SM
       SM=ZERO
       DO 36 M = 0,N-1
   36  SM=SM+V(M)*C(N-M,P)
       WGPLG=SGN(N)*SK+SGN(P)*(SM+V(N+P))
       RETURN
      END IF

      IF(X .GT. HALF) THEN
       X1=SGN(0)-X
       H=C1*X1+C2
       ALFA=H+H
       V(0)=SGN(0)
       U(0)=SGN(0)
       V(1)=LOG(DCMPLX(X1,ZERO))
       U(1)=LOG(X)
       DO 23 L = 2,P
   23  V(L)=V(1)*V(L-1)/L
       DO 26 L = 2,N
   26  U(L)=U(1)*U(L-1)/L
       SK=ZERO
       DO 24 K = 0,N-1
       P1=N-K
       R=X1**P1/FCT(P1)
       SM=ZERO
       DO 25 M = 0,P-1
       N1=P-M
       L=INDEX(10*N1+P1-10)
       B1=ZERO
       B2=ZERO
       DO 12 I = NC(L),0,-1
       B0=A(I,L)+ALFA*B1-B2
       B2=B1
   12  B1=B0
       Q=SGN(M)*(B0-H*B2)*R/P1**N1
   25  SM=SM+V(M)*Q
   24  SK=SK+U(K)*(S1(P1,P)-SM)
       WGPLG=SK+SGN(P)*U(N)*V(P)
       RETURN
      END IF

      L=INDEX(10*N+P-10)
      H=C1*X+C2
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO 11 I = NC(L),0,-1
      B0=A(I,L)+ALFA*B1-B2
      B2=B1
   11 B1=B0
      WGPLG=(B0-H*B2)*X**P/(FCT(P)*P**N)
      RETURN
 1000 FORMAT(/' ***** CERN SUBROUTINE WGPLG ... ILLEGAL VALUES',
     1        '   N = ',I3,'   P = ',I3)
      END

C-}}}



c..   
c..   Coefficient of the delta-function at NNLO.
c..   

      double precision function delta2(lfr,lfh)
      implicit none

      double precision lfh,lfr,lfac,lfren
      double precision z2,z3,pi
      double precision consCF,consCA,consNF
      double precision cf,ca,nf

      cf = 4.d0/3.d0
      ca = 3.d0
      nf = 5.d0
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0

      lfren = -lfh+lfr
      lfac = -lfh
      
      pi = 3.1415926535897932385d0
 
      consCF = 8.0d0/5.0d0 * z2**2 - 60.0d0 * z3 - 70.0d0 * z2 +
     .     511.0d0/4.0d0
      consCF = consCF + (18.0d0 - 32.0d0 * z2 ) * lfac**2 +
     .       ( -93.0d0 + 24.0d0 * z2 + 176.0d0 * z3 ) * lfac

      consCA = -1535.0d0/12.0d0 + 592.0d0/9.0d0 * z2 + 28.0d0 * z3
     .     -12.0d0/5.0d0 * z2**2
      consCA = consCA + 11.0d0*lfac**2 - 22.0d0*lfac*lfren +
     .       ( 17.0d0/3.0d0  + 88.0d0/3.0d0 * z2 - 24.0d0 * z3)*lfac +
     .       ( 176.0d0/3.0d0 - 88.0d0/3.0d0 * z2 )*lfren

      consNF = 8.0d0 * z3 - 112.0d0/9.0d0 * z2 + 127.0d0/6.0d0
      consNF = consNF - 2.0d0*lfac**2 + 4.0D0*lfac*lfren
     .       - (  2.0d0/3.0d0 + 16.0d0/3.0d0 * z2 ) * lfac +
     .         (-32.0d0/3.0d0 + 16.0d0/3.0d0 * z2 ) * lfren
      delta2 = cf**2 * consCF + cf * ca * consCA + nf*cf * consNF
      delta2 = delta2/(4.d0*pi)**2
      return
      end

c..
c..   The plus-distributions at NNLO.
c..   
      double precision function dterms2(xt,lfh,lfr)
      implicit none
      double precision cf,nf,ca,pi
      double precision partCF,partCA,partNF
      double precision xt,lfh,lfr,lfac,lfren
      double precision log1mz,dd0,dd1,dd2,dd3
      double precision z2,z3,z4,z
      z=xt

      cf = 4.d0/3.d0
      nf = 5.d0
      ca = 3.d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0

      lfren = -lfh+lfr
      lfac = -lfh      
      pi = 3.1415926535897932385d0

      log1mz = dlog(1.0d0-z)
      dd0 = 1.0d0/(1.0d0-z)
      dd1 = dd0 * log1mz
      dd2 = dd0 * log1mz**2
      dd3 = dd0 * log1mz**3

      partCF = 128.0d0 * dd3 - (128.0d0*z2 + 256.0d0 )*dd1 +
     .     256.0d0 * z3 * dd0
      partCF = partCF + 192.0d0 * lfac * dd2 +
     .       ( 64.0d0 * lfac + 96.0d0 ) * lfac * dd1 +
     .       ( 48.0D0 * lfac - 128.0d0 - 64.0d0*z2 ) * lfac * dd0

      partCF = cf**2 * partCF

c--- CA
      partCA = -176.0d0/3.0d0 * dd2 + ( 1072.0d0/9.0d0 - 32.0d0*z2)*dd1+
     .     ( 56.0d0*z3 + 176.0d0/3.0d0*z2 -1616.0d0/27.0d0 )*dd0
      partCA = partCA - 176.0d0/3.0d0*lfren*dd1 +
     .       ( 44.0d0/3.0d0*lfac**2 - 88.0d0/3.0d0*lfac*lfren +
     .         ( 536.0d0/9.0d0 - 16.0d0*z2 )*lfac ) * dd0

      partCA = ca*cf * partCA

c -- NF

      partNF = 32.0d0/3.0d0 * dd2 - 160.0d0/9.0d0 * dd1 +
     .     ( 224.0d0/27.0d0 - 32.0d0/3.0d0*z2 )*dd0
      partNF = partNF + 32.0d0/3.0d0*lfren * dd1 +
     .       ( -8.0d0/3.0d0*lfac**2 + 16.0d0/3.0d0*lfac*lfren
     .         -80.0d0/9.0d0*lfac ) * dd0
      partNF = nf*cf * partNF

      dterms2 = partCF  + partCA + partNF
      dterms2 = dterms2/(4.d0*pi)**2
      return
      end

c.......................................................................
c............  Unchanged from bbh .......................................
c.......................................................................
      double precision function dtsub2(lfh,lfr,tauh)
c..
c..   Contributions arising from the fact that the integrals over
c..   plus-distributions do not run from 0 to 1, but from z to 1.
c..   
      implicit double precision (a-z)
      dtsub2 = 0.
      end



C........................Running mass (Marius) ....................
      double precision function runmass(mass0,api0,apif,nf,nloop)
c..
c..   evaluates the running of the MS-bar quark mass
c..   by expanding the equation
c..   
c..   m(mu) = m(mu0) * exp( \int_a0^af dx gammam(x)/x/beta(x) )
c..   
c..   in terms of alpha_s. The results agree with RunDec.m.
c..   
c..   
c..   Input:
c..   ------
c..   mass0  :  m(mu0)
c..   api0   :  alpha_s(mu0)/pi
c..   apif   :  alpha_s(muf)/pi
c..   nf     :  number of flavors
c..   nloop  :  order of calculation (nloop=1..4)
c..
c..   Output:
c..   -------
c..   massout:  m(muf)
c..   
      implicit real*8 (a-h,o-z)
      real*8 mass0,massfun
      real*8 nf
      external massfun
      parameter(accmass=1.d-6)
      data z3/1.2020569031595942853997/,
     &     z5/1.0369277551433699263/,
     &     pi/3.1415926535897932381/

      beta0 = (33 - 2*nf)/12.d0
      beta1 = (102 - (38*nf)/3.d0)/16.d0
      beta2 = (2857/2.d0 - (5033*nf)/18.d0 + (325*nf**2)/54.d0)/64.d0
      beta3 = (149753/6.d0 + (1093*nf**3)/729.d0 + 3564*z3 + nf**2
     &     *(50065/162.d0 + (6472*z3)/81.d0) - nf*(1078361/162.d0 +
     &     (6508*z3)/27.d0))/256.d0
      
      gamma0 = 1.d0
      gamma1 = (67.33333333333333d0 - (20*nf)/9.d0)/16.d0
      gamma2 = (1249.d0 - (140*nf**2)/81.d0 + 2*nf*(-20.59259259259259d0
     &     - 48*z3) +(8*nf*(-46 + 48*z3))/9.d0)/64.d0
      gamma3 = (28413.91975308642d0 + (135680*z3)/27.d0 + nf**3*(-1
     &     .3662551440329218d0 + (64*z3)/27.d0) + nf**2*(21
     &     .57201646090535d0 - (16*Pi**4)/27.d0 + (800*z3)/9.d0) - 8800
     &     *z5 + nf*(-3397.1481481481483d0 + (88*Pi**4)/9.d0 - (34192
     &     *z3)/9.d0 + (18400*z5)/9.d0))/256.d0
      

      bb1 = beta1/beta0
      bb2 = beta2/beta0
      bb3 = beta3/beta0

      cc0 = gamma0/beta0
      cc1 = gamma1/beta0
      cc2 = gamma2/beta0
      cc3 = gamma3/beta0

      cfunc1 = 1.d0
      cfunc2 = cc1 - bb1*cc0
      cfunc3 = 1/2.d0*((cc1-bb1*cc0)**2 + cc2 - bb1*cc1 + bb1**2*cc0 -
     &     bb2*cc0)
      cfunc4 = (1/6*(cc1 - bb1*cc0)**3 + 1/2*(cc1 - bb1*cc0)*(cc2 - bb1
     &     *cc1 + bb1**2*cc0 - bb2*cc0) + 1/3*(cc3 - bb1*cc2 + bb1**2
     &     *cc1 - bb2*cc1 - bb1**3*cc0 + 2*bb1*bb2*cc0 - bb3*cc0))

      if (nloop.lt.4) then
         cfunc4 = 0.d0
         if (nloop.lt.3) then
            cfunc3 = 0.d0
            if (nloop.lt.2) then
               cfunc2 = 0.d0
               if (nloop.lt.1) then
                  cfunc1 = 0.d0
               endif
            endif
         endif
      endif

      cfuncmu0 = cfunc1 + cfunc2*api0 + cfunc3*api0**2 + cfunc4*api0**3
      cfuncmuf = cfunc1 + cfunc2*apif + cfunc3*apif**2 + cfunc4*apif**3

      runmass = mass0*(apif/api0)**cc0*cfuncmuf/cfuncmu0
      end


      DOUBLE PRECISION FUNCTION CQBCAR(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr,z,help,zmin
      double precision logz,log1mz,di1mz,s1mz,dimz,trimz,smz
      double precision log1pz,tri1mz
      double precision trimco,tripco
      double precision part1,part2,scale
      double precision z2,z3,z4
      double precision ca,cf,nf
      double precision dilog,trilog,wgplg
      external dilog,trilog,wgplg
      
      ca=3.d0
      cf=4.d0/ca
      nf=5d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      Z = xx

C.
      LOGZ=DLOG(Z)
      ZMIN=1.0D0-Z
      LOG1MZ=DLOG(ZMIN)
      DI1MZ=WGPLG(1,1,ZMIN)
      TRI1MZ=WGPLG(2,1,ZMIN)
      S1MZ=WGPLG(1,2,ZMIN)
      HELP=( 1.0D0 + Z ) *
     .          ( 20.0D0*S1MZ - 28.0D0*z3 +
     .            16.0D0*LOG1MZ*DI1MZ - 8.0D0*LOGZ*DI1MZ +
     .            16.0D0*LOG1MZ*z2 - 8.0D0*LOGZ*z2 +
     .            88.0D0/3.0D0*LOG1MZ**2 ) +
     .     (-8.0D0*S1MZ - 24.0D0*TRI1MZ
     .      -16.0D0*LOG1MZ*DI1MZ + 16.0D0*LOGZ*DI1MZ +
     .       16.0D0*LOGZ*z2 + 8.0D0/3.0D0*DI1MZ +
     .       280.0D0/3.0D0*LOG1MZ*LOGZ - 29.0D0*LOGZ**2
     .      -208.0D0/3.0D0*LOGZ )/ZMIN
     .      -DI1MZ*( 32.0D0/3.0D0 + 56.0D0/3.0D0*Z )
     .      -z2*( 76.0D0/3.0D0 + 100.0D0/3.0D0*Z )
     .      -LOG1MZ*LOGZ*( 176.0D0/3.0D0 + 176.0D0/3.0D0*Z ) +
     .       LOGZ**2*( 55.0D0/3.0D0 + 55.0D0/3.0D0*Z )
     .      -LOG1MZ*( 152.0D0/9.0D0 + 956.0D0/9.0D0*Z ) +
     .       LOGZ*( 52.0D0/3.0D0 + 218.0D0/3.0D0*Z )
     .      -446.0D0/27.0D0 + 2278.0D0/27.0D0*Z

      HELP = HELP +
     .     lfh**2*(  - 22.D0/3.D0 - 22.D0/3.D0*Z ) +
     .     lfh*lfr * ( 44.D0/3.D0 + 44.D0/3.D0*Z ) +
     .     lfh   *( LOGZ * (  - 44.D0/3.D0 - 44.D0/3.D0*Z +
     .     52.D0/3.D0/ZMIN )
     .     + DI1MZ * ( 8.D0 + 8.D0*Z - 16.D0/ZMIN )
     .     + z2 * ( 8.D0 + 8.D0*Z )
     .     - 76/9.D0 - 496.D0/9.D0*Z                   ) +
     .     lfr   *( LOGZ * (  - 44.D0/3.D0 - 44.D0/3.D0*Z +
     .     88.D0/3.D0/ZMIN )
     .     + LOG1MZ * ( 88.D0/3.D0 + 88.D0/3.D0*Z )    )
      CQBCAR=CA*CF*HELP
      RETURN
      END

      DOUBLE PRECISION FUNCTION CQBNF(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr,z,help
      double precision logz,log1mz,s1mz,dimz,trimz,smz
      double precision d0,d1,d2
      double precision log1pz,di1mz
      double precision trimco,tripco
      double precision part1,part2,scale
      double precision z2,z3,z4
      double precision ca,cf,nf
      double precision dilog,trilog,wgplg
      external dilog,trilog,wgplg
      
      ca=3.d0
      cf=4.d0/ca
      nf=5d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      Z = xx
      LOG1MZ=DLOG(1.0D0-Z)
      D0=1.0D0/(1.0D0-Z)
      D1=D0*LOG1MZ
      D2=D1*LOG1MZ
      HELP=32.0D0/3.0D0*D2 - 160.0D0/9.0D0*D1 +
     .     ( 224.0D0/27.0D0 - 32.0D0/3.0D0*z2 )*D0
      HELP=HELP + 32.0D0/3.0D0*lfr*D1 +
     .     ( -8.0D0/3.0D0*lfh**2 + 16.0D0/3.0D0*lfh*lfr
     .     -80.0D0/9.0D0*lfh )*D0
     
      CQBNF=NF*CF * HELP
      RETURN
      END

      DOUBLE PRECISION FUNCTION CQBNFR(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr,z,zmin,help
      double precision logz,s1mz,dimz,trimz,smz,log1mz
      double precision log1pz,di1mz
      double precision trimco,tripco
      double precision part1,part2,scale
      double precision z2,z3,z4
      double precision ca,cf,nf
      double precision dilog,trilog,wgplg
      external dilog,trilog,wgplg
      
      ca=3.d0
      cf=4.d0/ca
      nf=5d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      Z = xx
      LOGZ=DLOG(Z)
      ZMIN=1.0D0-Z
      LOG1MZ=DLOG(ZMIN)
      DI1MZ=WGPLG(1,1,ZMIN)
      HELP=( 1.0D0 + Z ) * ( 16.0D0/3.0D0*z2 + 8.0D0/3.0D0*DI1MZ -
     .             10.0D0/3.0D0*LOGZ**2 + 32.0D0/3.0D0* LOG1MZ*LOGZ -
     .             16.0D0/3.0D0*LOG1MZ**2 ) +
     .     ( -64.0D0/3.0D0*LOG1MZ*LOGZ - 8.0D0/3.0D0*DI1MZ +
     .        8.0D0*LOGZ**2 + 40.0D0/3.0D0*LOGZ ) / ZMIN +
     .     LOG1MZ * ( - 16.0D0/9.0D0 + 176.0D0/9.0D0*Z ) +
     .     LOGZ * ( - 4.0D0/3.0D0 - 44.0D0/3.0D0*Z ) +
     .     188.0D0/27.0D0 - 412.0D0/27.0D0*Z
      HELP = HELP +
     .     lfh**2*( 4.D0/3.D0 + 4.D0/3.D0*Z ) +
     .     lfh*lfr * (  - 8.D0/3.D0 - 8.D0/3.D0*Z ) +
     .     lfh   *( LOGZ * ( 8.D0/3.D0 + 8.D0/3.D0*Z
     .     - 16.D0/3.D0/ZMIN )
     .     - 8.D0/9.D0 + 88.D0/9.D0*Z                  ) +
     .     lfr   *( LOGZ * ( 8.D0/3.D0 + 8.D0/3.D0*Z
     .     - 16.D0/3.D0/ZMIN )
     .     + LOG1MZ * (  - 16.D0/3.D0 - 16.D0/3.D0*Z ) )

      CQBNFR=NF*CF*HELP
      RETURN
      END


      DOUBLE PRECISION FUNCTION CQBCFR(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr,z,zmin,help
      double precision logz,s1mz,dimz,trimz,smz,log1mz
      double precision log1pz,di1mz,tri1mz
      double precision trimco,tripco
      double precision part1,part2,scale
      double precision z2,z3,z4
      double precision ca,cf,nf
      double precision dilog,trilog,wgplg
      external dilog,trilog,wgplg
      
      ca=3.d0
      cf=4.d0/ca
      nf=5d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      Z = xx
      LOGZ=DLOG(Z)
      ZMIN=1.0D0-Z
      LOG1MZ=DLOG(ZMIN)
      DI1MZ=WGPLG(1,1,ZMIN)
      TRI1MZ=WGPLG(2,1,ZMIN)
      S1MZ=WGPLG(1,2,ZMIN)
      HELP=( 1.0D0 + Z ) *
     .          ( 48.0D0*S1MZ - 32.0D0*TRI1MZ - 128.0D0*z3 +
     .            24.0D0*LOG1MZ*DI1MZ + 24.0D0*LOGZ*DI1MZ +
     .            64.0D0*LOG1MZ*z2 -96.0D0*LOGZ*z2
     .           -64.0D0*LOG1MZ**3 + 156.0D0*LOG1MZ**2*LOGZ
     .           -96.0D0*LOG1MZ*LOGZ**2 + 50.0D0/3.0D0*LOGZ**3 ) +
     .     (-64.0D0*S1MZ - 16.0D0*TRI1MZ +
     .       48.0D0*LOG1MZ*DI1MZ - 48.0D0*LOGZ*DI1MZ +
     .       128.0D0*LOGZ*z2
     .      -248.0D0*LOG1MZ**2*LOGZ + 144.0D0*LOG1MZ*LOGZ**2
     .      -24.0D0*LOGZ**3 + 112.0D0*LOGZ )/ZMIN +
     .     DI1MZ*( -24.0D0 - 32.0D0*Z ) + z2*( 64.0D0 - 64.0D0*Z ) +
     .     LOG1MZ**2*( - 64.0D0 + 64.0D0*Z ) +
     .     LOG1MZ*LOGZ*( 64.0D0 - 112.0D0*Z ) +
     .     LOGZ**2*( - 8.0D0 + 24.0D0*Z ) +
     .     LOG1MZ*( 256.0D0 + 12.0D0*Z ) +
     .     LOGZ*( - 104.0D0 + 48.0D0*Z ) - 72.0D0 + 48.0D0*Z

      HELP = HELP +
     .     LFH**2*(  LOGZ * ( 24.D0 + 24.D0*Z - 32.D0/ZMIN )
     .     + LOG1MZ * (  - 32.D0 - 32.D0*Z )
     .     - 40.D0 - 8.D0*Z                          ) +
     .     LFH   *( LOGZ**2 * (  - 36.D0 - 36.D0*Z + 48.D0/ZMIN )
     .     + LOGZ*LOG1MZ * ( 144.D0 + 144.D0*Z - 224.D0/ZMIN )
     .     + LOGZ * ( 56.D0 - 24.D0*Z - 48.D0/ZMIN )
     .     + LOG1MZ**2 * (  - 96.D0 - 96.D0*Z )
     .     + LOG1MZ * (  - 112.D0 + 16.D0*Z )
     .     + DI1MZ * ( 16.D0 + 16.D0*Z + 32.D0/ZMIN )
     .     + z2 * ( 32.D0 + 32.D0*Z )
     .     + 120.D0 + 16.D0*Z                                )

      CQBCFR=CF*CF*HELP
      RETURN
      END



      DOUBLE PRECISION FUNCTION CQBACD(xx,lfh,lfr)
      implicit none
      double precision xx,lfh,lfr,z,zmin,help
      double precision logz,s1mz,dimz,trimz,smz,log1mz
      double precision log1pz,di1mz,tri1mz
      double precision trimco,tripco
      double precision part1,part2,scale
      double precision z2,z3,z4
      double precision ca,cf,nf
      double precision delacd,frterm
      double precision dilog,trilog,wgplg
      external dilog,trilog,wgplg
      
      ca=3.d0
      cf=4.d0/ca
      nf=5d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      Z = xx
      LOGZ=DLOG(Z)
      ZMIN=1.0D0-Z
      LOG1MZ=DLOG(ZMIN)
      DI1MZ=WGPLG(1,1,ZMIN)
      TRI1MZ=WGPLG(2,1,ZMIN)
      S1MZ=WGPLG(1,2,ZMIN)
      DELACD=CF * ( CF - CA/2 ) * (
     . ( LOGZ**3 * 16.D0/3.D0  - LOGZ**2*LOG1MZ * 16.D0
     . + LOGZ**2 * 15.D0       - LOGZ*LOG1MZ * 24.D0
     . - LOGZ*DI1MZ * 24.D0   + LOGZ * 24.D0
     . - LOG1MZ*DI1MZ * 32.D0 - DI1MZ * 12.D0
     . + TRI1MZ * 32.D0       - S1MZ * 72.D0          )/ZMIN
     . + LOGZ**3 * (  - 2.D0 - 2.D0*Z )
     . + LOGZ**2*LOG1MZ * ( 8.D0 + 8.D0*Z )
     . + LOGZ**2 * ( 4.D0 + 4.D0*Z )
     . + LOGZ*LOG1MZ * (  - 16.D0 - 16.D0*Z )
     . + LOGZ*DI1MZ * ( 16.D0 + 16.D0*Z )
     . + LOGZ * ( 32.D0 - 30.D0*Z )
     . + LOG1MZ*DI1MZ * ( 16.D0 + 16.D0*Z )
     . + LOG1MZ * (  - 64.D0 + 56.D0*Z )
     . + DI1MZ * (  - 20.D0 - 20.D0*Z )
     . + TRI1MZ * (  - 24.D0 - 24.D0*Z )
     . + S1MZ * ( 36.D0 + 36.D0*Z )
     . + 94.D0 - 78.D0*Z )
      FRTERM=CF * ( CF - CA/2 ) * lfh * (
     .     LOGZ**2 * ( 4.D0 + 4.D0*Z - 8.D0/ZMIN )
     .     + LOGZ * (  - 8.D0 - 8.D0*Z - 12.D0/ZMIN )
     .     + DI1MZ * ( 8.D0 + 8.D0*Z - 16.D0/ZMIN )
     .     - 32.D0 + 28.D0*Z                 )

      CQBACD=DELACD+FRTERM
      RETURN
      END

      double precision function CQQCRI(xx,lfh,lfr)
c--- CQQCRI
      implicit none
      double precision partCRI
      double precision xx,lfh,lfr,z
      double precision xm1,xp1,dlnx,dlxm1,dlxp1
      double precision dli2a,dli3a,zcom,s1mz,dimz,trimz,smz
      double precision trimco,tripco
      double precision part1,part2,scale
      double precision z2,z3,z4
      double precision ca,cf,nf
      double precision dilog,trilog,wgplg
      double precision cax,cv,cd,sw
      external dilog,trilog,wgplg

      sw  = 0.2228972225239183d0
      cax = 1.d0
      cv  = -1.d0 + 4.d0/3.d0*sw
      cd  = 1.d0 + cv**2

      
      ca=3.d0
      cf=4.d0/ca
      nf=5d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)
      xp1   = 1.0d0 + xx
      dlxp1 = dlog(xp1)
      dli2a = dilog(xm1)
      dli3a = trilog(xm1)
      dimz  = wgplg(1,1,-xx)
      trimz = wgplg(2,1,-xx)
      zcom  = xm1/xp1
      tripco= wgplg(2,1,zcom)
      trimco= wgplg(2,1,-zcom)
      s1mz  = wgplg(1,2,xm1)
      smz   = wgplg(1,2,-xx)

      Z = xx

      partCRI = CF * ( CF - CA/2 ) * (
     .   dlnx**3 * (  - 4.D0/3.D0 - 4.D0/3.D0*Z**2 + 8.D0/3.D0*Z )
     . + dlnx**2 * (  - 6.D0 - 6.D0*Z**2 + 12.D0*Z )
     . + dlnx*dli2a * (  - 8.D0 - 8.D0*Z**2 + 16.D0*Z )
     . + dlnx * (  - 14.D0 + 12.D0*Z )
     . + dli2a * (  - 12.D0 - 12.D0*Z**2 + 24.D0*Z )
     . + dli3a * ( 8.D0 + 8.D0*Z**2 - 16.D0*Z )
     . + S1MZ * (  - 8.D0 - 8.D0*Z**2 + 16.D0*Z )
     .     - 15.D0 - 13.D0*Z**2 + 28.D0*Z )
      CQQCRI = partCRI
      return
      end

      double precision function CQQCRF(xx,lfh,lfr)
      implicit none
      double precision partCRF
      double precision xx,lfh,lfr,z
      double precision xm1,xp1,dlnx,dlxm1,dlxp1
      double precision dli2a,dli3a,zcom,s1mz,dimz,trimz,smz
      double precision trimco,tripco
      double precision part1,part2,scale
      double precision z2,z3,z4
      double precision ca,cf,nf
      double precision dilog,trilog,wgplg
      double precision cax,cv,cd,sw
      external dilog,trilog,wgplg

      sw  = 0.2228972225239183d0
      cax = 1.d0
      cv  = -1.d0 + 4.d0/3.d0*sw
      cd  = 1.d0 + cv**2

      
      ca=3.d0
      cf=4.d0/ca
      nf=5d0
      
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0
      z4 = 1.0823232337111381915d0
      
      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)
      xp1   = 1.0d0 + xx
      dlxp1 = dlog(xp1)
      dli2a = dilog(xm1)
      dli3a = trilog(xm1)
      dimz  = wgplg(1,1,-xx)
      trimz = wgplg(2,1,-xx)
      zcom  = xm1/xp1
      tripco= wgplg(2,1,zcom)
      trimco= wgplg(2,1,-zcom)
      s1mz  = wgplg(1,2,xm1)
      smz   = wgplg(1,2,-xx)

      Z = xx


      partCRF =
     .     dlnx**3 * ( 2.D0 - 2.D0*Z - 16.D0/3.D0/xp1 )
     .     + dlnx**2*dlxm1 * (  - 8.D0 + 8.D0*Z + 16.D0/xp1 )
     .     + dlnx**2*dlxp1 * (  - 24.D0 + 24.D0*Z + 56.D0/xp1 )
     .     + dlnx**2 * (  - 4.D0 - 12.D0*Z )
     .     + dlnx*dlxm1*dlxp1 * ( 32.D0 - 32.D0*Z - 64.D0/xp1 )
     .     + dlnx*dlxm1 * ( 16.D0 + 16.D0*Z )
     .     + dlnx*dlxp1**2 * (  - 16.D0/xp1 )
     .     + dlnx*dlxp1 * ( 8.D0 + 8.D0*Z )
     .     + dlnx*dli2a * (  - 24.D0 + 24.D0*Z + 48.D0/xp1 )
     .     + dlnx*DIMZ * (  - 32.D0 + 32.D0*Z + 64.D0/xp1 )
     .     + dlnx*z2 * (  - 8.D0 + 8.D0*Z + 24.D0/xp1 )
     .     + dlnx * (  - 18.D0 + 14.D0*Z )
     .     + dlxm1*DIMZ * ( 32.D0 - 32.D0*Z - 64.D0/xp1 )
     .     + dlxm1*z2 * ( 16.D0 - 16.D0*Z - 32.D0/xp1 )
     .     + dlxm1 * ( 32.D0 - 32.D0*Z )
     .     + dlxp1*DIMZ * (  - 32.D0/xp1 )
     .     + dlxp1*z2 * (  - 16.D0/xp1 )
      partCRF = partCRF
     .     + dli2a * ( 24.D0 + 8.D0*Z )
     .     + DIMZ * ( 8.D0 + 8.D0*Z )
     .     + z2 * ( 4.D0 + 4.D0*Z )
     .     + dli3a * ( 32.D0 - 32.D0*Z - 64.D0/xp1 )
     .     + TRIMZ * ( 16.D0 - 16.D0*Z - 16.D0/xp1 )
     .     + TRIMCO * ( 32.D0 - 32.D0*Z - 64.D0/xp1 )
     .     + TRIPCO * (  - 32.D0 + 32.D0*Z + 64.D0/xp1 )
     .     + S1MZ * (  - 32.D0 + 32.D0*Z + 64.D0/xp1 )
     .     + SMZ * (  - 32.D0/xp1 )
     .     + z3 * ( 12.D0 - 12.D0*Z - 8.D0/xp1 )
     .     - 34.D0 + 34.D0*Z

      partCRF = partCRF + lfh * (
     .     dlnx**2 * (  - 4.D0 + 4.D0*Z + 8.D0/xp1 )
     .     + dlnx*dlxp1 * ( 16.D0 - 16.D0*Z - 32.D0/xp1 )
     .     + dlnx * ( 8.D0 + 8.D0*Z )
     .     + DIMZ * ( 16.D0 - 16.D0*Z - 32.D0/xp1 )
     .     + z2 * ( 8.D0 - 8.D0*Z - 16.D0/xp1 )
     .     + 16.D0 - 16.D0*Z                 )      
      partCRF = cf * ( cf - ca/2 ) * partCRF
      CQQCRF = partCRF
      return
      end
