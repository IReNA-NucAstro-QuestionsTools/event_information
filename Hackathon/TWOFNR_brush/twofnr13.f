******************************************************************
*
*     twofnr5.f includes the changes required for the data
*     sets produced by front5z.f, as above. This version also
*     calculates the 'ANC' [actually the square of the ANC]
*     of the single particle orbital (assuming at present 
*     that the transferred particle is neutral, i.e. (d,p), 
*     (p,d) etc.). It also prints the rms radius, as before. 
*     Screen output has been minimised. 
*
*     Version twofnr6.f includes the mods (started December 05)
*     to rotate the transfer amplitudes to alternative coordinate 
*     systems for m-substate populations/cross sections - WNC.  
*
*     Version 10 unchanged but renumbered to keep pace with the
*     upwardly mobile front end numbering - now at 10.
*
*     This code has fundamental constants adjusted to agree with 
*     those used in FRESCO for better checking capabilities etc.
*     This version also uses the coulfg Coulomb functions routine
*     for improved accuracy and to help track down an apparent 
*     instability in calculations at higher incident energies. 
*     (NSCL January 2012).
*
*     Changed amplitude storage arrangements for reading/mixing
*     Redimensioned arrays for larger final state spins, WNC 2012
*
c     *****( twostp main )*****                                         
      implicit real*8(a-h,o-z)
      character*80 infname                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/pot/vd(2),wd(2),vsod(2),wsod(2),rrd(2),ard(2),rcd(2),      
     1csdgd(2),rid(2),aid(2),rsord(2),asord(2),rsoid(2),asoid(2),       
     2ktop(2),rgd(2),agd(2),ktisp(2),risrd(2),aisrd(2),visd(2),risid(2),
     3aisid(2),wisd(2),csd,drhod(12),rhobn,rhoin,rhosor,rhosoi,rhoisr,  
     4rhoisi,rhobc,rhobng,fkayar,fkayai,fkasor,fkasoi,fkaisr,fkaisi,    
     5fkaybg,vbe,wbe,vsoe,wsoe,vise,wise,vbcd(905,2),wbcd(905,2),       
     6vsov(905,2),wsov(905,2),pnloc(2)                                  
      common/coul/z(16),sigld(90,2),fd(90,2),gd(90,2),fpd(90,2),        
     1gpd(90,2)                                                         
      common/theta/jmax,j1min,j1max,ktheta,mmax,nnnn,
     1thetad(181),sint(181),    
     1cost(181),plmd(181,10),xdwb(181),poldw(181),asymdw(181),
     2tenpol(181,6),plmq(181,10)
      common/ffc/amp,d02,fts,ta,tb,fnrng,ktff(5),nrange,nafac,nbfac,lsa 
     1,lsb,iread(15),jjj,fread(28),wread(21),ktrl(20),ff(905),ffi(905), 
     2ffone(905),nrangd(2,4),n1xd(2,4),itap(66,10,2),frfac(2,4),frso(2) 
      complex*16 ffone,frfac,samp,samf                                  
      common/clebma/faclog(500)                                         
      common/add/itape                                                  
      common/zf/ktzf(2),kk2,kk3                                         
      common/out/idknt                                                  
      common/bamp/samp(181,10,3,3),samf(181,10,3,3)     
      print*,'input file name ?'
      read '(a)', infname
      open(35,file=infname,status='old')
      iskip1=0                                                          
      iskip2=0                                                          
      mxsi=0                                                            
      mxsf=0                                                            
      mxjp=0                                                            
      do 999 ijt=1,2                                                    
      ktzf(ijt)=0                                                       
      faclog(ijt)=0.d0                                                  
  999 beta(ijt)=0.0d0                                                   
      kk2=0                                                             
      kk3=0                                                             
      do 1000 ijt=1,5                                                   
 1000 numrun(ijt)=0                                                     
      mmax=0                                                            
      call null                                                         
      fn=1.0d0                                                          
      do 200 i=3,500                                                    
      fn=fn+1.0d0                                                       
  200 faclog(i)=faclog(i-1)+log(fn)                                     
      noline=55                                                         
      itape=0                                                           
 2000 npgs=0                                                            
      kcheck=0                                                          
      iskip1=0                                                          
      numrun(5)=numrun(5)+1                                             
      do 100 n1=1,905                                                   
  100 ffone(n1)=(0.,0.)                                                 
      do 120 n1=1,181                                                   
      xdwb(n1)=0.0                                                      
      poldw(n1)=0.0                                                     
  120 asymdw(n1)=0.0                                                    
      call basic                                                        
  805 call ffgenz(1,2,2,1)                                              
      do 810 nr=1,nrmax                                                 
  810 ffone(nr)=cmplx(ff(nr),ffi(nr))+ffone(nr)                         
      if(ktff(4).gt.0) go to 805                                        
      call tasktr                                                       
      go to 2000                                                        
      end                                                               
c     ******************************************************************
      subroutine basic                                                  
      implicit real*8(a-h,o-z)                                          
*     dimension aaa(90),bbb(90)
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/pot/vd(2),wd(2),vsod(2),wsod(2),rrd(2),ard(2),rcd(2),      
     1csdgd(2),rid(2),aid(2),rsord(2),asord(2),rsoid(2),asoid(2),       
     2ktop(2),rgd(2),agd(2),ktisp(2),risrd(2),aisrd(2),visd(2),risid(2),
     3aisid(2),wisd(2),csd,drhod(12),rhobn,rhoin,rhosor,rhosoi,rhoisr,  
     4rhoisi,rhobc,rhobng,fkayar,fkayai,fkasor,fkasoi,fkaisr,fkaisi,    
     5fkaybg,vbe,wbe,vsoe,wsoe,vise,wise,vbcd(905,2),wbcd(905,2),       
     6vsov(905,2),wsov(905,2),pnloc(2)                                  
      common/coul/z(16),sigld(90,2),fd(90,2),gd(90,2),fpd(90,2),        
     1gpd(90,2)                                                         
      common/jat/ sigmad(90),f(90),fp(90),g(90),gp(90)                  
      call inputa                                                      
      fnrmax=nrmax                                                      
      dr=rmax/fnrmax                                                    
      ichmx=2+nubchn                                                    
      do  500  i=1,ichmx                                                
      ichnl=i                                                           
      csd=csdgd(i)                                                      
      tautau(i)=(tmas(i)-2.0*tz(i))*(pmas(i)-2.0*pz(i))*0.25/(tmas(i)   
     1*pmas(i))                                                         
      t1=pmas(i)+tmas(i)                                                
      fmu=pmas(i)*tmas(i)/t1                                            
      fmud(i)=fmu                                                       
      if(i-1) 1,1,2                                                     
    1 ecm=elabi*tmas(1)/t1                                              
      ecmi=ecm                                                          
      go  to  5                                                         
    2 ecm=ecmi+qvlue(i)                                                 
    5 ecmd(i)=ecm                                                       
      ecm=abs(ecm)                                                      
*     fkay=0.2195376*sqrt(fmu*ecm)  
      fkay=0.218735d0*sqrt(fmu*ecm)                                     
      fk(i)=fkay                                                        
*     eta=0.15805086*zz(i)*sqrt(fmu/ecm)    
      eta=0.1574855d0*zz(i)*sqrt(fmu/ecm)                             
      drd(i)=dr*tmas(1)/tmas(i)                                         
      drho=fkay*drd(i)                                                  
      drhod(i)=drho                                                     
      drho56=0.83333333*(drho*drho)                                     
      drhod(i+8)=drho56                                                 
      drhod(i+4)=drho56*0.1                                             
      t2=fkay*tmas(i)**0.3333333                                        
      t3=fkay*(pmas(i)**0.33333333+tmas(i)**0.33333333)                 
      rhobn=rrd(i)*t2                                                   
      if(rhobn.lt.0.) rhobn=-rrd(i)*t3                                  
      rhoin=rid(i)*t2                                                   
      if(rhoin.lt.0.) rhoin=-rid(i)*t3                                  
      rhosor=rsord(i)*t2                                                
      if(rhosor.lt.0.) rhosor=-rsord(i)*t3                              
      rhosoi=rsoid(i)*t2                                                
      if(rhosoi.lt.0.) rhosoi=-rsoid(i)*t3                              
      rhoisr=risrd(i)*t2                                                
      if(rhoisr.lt.0.) rhoisr=-risrd(i)*t3                              
      rhoisi=risid(i)*t2                                                
      if(rhoisi.lt.0.) rhoisi=-risid(i)*t3                              
      rhobc=rcd(i)*t2                                                   
      if(rhobc.lt.0.) rhobc=-rcd(i)*t3                                  
      rhobng=rgd(i)*t2                                                  
      if(rhobng.lt.0.) rhobng=-rgd(i)*t3                                
      rhomax=drho*fnrmax                                                
      fkayar=ard(i)*fkay                                                
      fkayai=aid(i)*fkay                                                
      fkasor=asord(i)*fkay                                              
      fkasoi=asoid(i)*fkay                                              
      fkaisr=aisrd(i)*fkay                                              
      fkaisi=aisid(i)*fkay                                              
      fkaybg=agd(i)*fkay                                                
      vbe=vd(i)/ecm                                                     
      wbe=wd(i)/ecm                                                     
      vsoe=vsod(i)/ecm                                                  
      wsoe=wsod(i)/ecm                                                  
      vise=visd(i)/ecm                                                  
      wise=wisd(i)/ecm                                                  
      z(i)=eta                                                          
      z(i+12)=eta/(rcd(i)**3*tmas(i))                                   
      lmax=lmaxd(i)+1                                                   
      zzz=zz(i)                                                         
      if(zzz)  10,10,20                                                 
   10 do 15 l=1,lmax                                                    
      sigld(l,i)=0.0                                                    
   15 continue                                                          
      go to 40                                                          
   20 call sigma(lmax,eta,sigma0,sigmad)                                
      do 35 l=1,lmax                                                    
      sigld(l,i)=sigmad(l)                                              
   35 continue                                                          
   40 continue                                                          
      mode1=1
      kfn=0
      zero=0.d0
      xlmax=lmax
      call coulfg(rhomax,eta,zero,xlmax,f,g,fp,gp,mode1,kfn,ifail)
*     call coulfn(eta,rhomax,drho,lmax)                                 
*     do l=1,lmax
*     print*,aaa(l),fp(l),abs(aaa(l)-fp(l))
*     print*,bbb(l),gp(l),abs(bbb(l)-gp(l))
*     enddo
      do 55 l=1,lmax                                                    
      fd(l,i)=f(l)                                                      
      gd(l,i)=g(l)                                                      
      fpd(l,i)=fp(l)                                                    
      gpd(l,i)=gp(l)                                                    
   55 continue                                                          
      nr3max=nrmax+3                                                    
      nr3min=nrmax-3                                                    
      call  pgen2                                                       
  500 continue                                                          
      call writea                                                       
      return                                                            
      end                                                               
c     ******************************************************************
      block data                                                        
      implicit real*8(a-h,o-z)                                          
      common/zrff/ra(15),ia,ib,rc(919)                                  
      data ra,ia,ib,rc/15*0.0,0,0,919*0.0/                              
      end                                                               
c     ******************************************************************
      function  cleb(ia,id,ib,ie,ic,if)                                 
      implicit real*8(a-h,o-z)                                          
      common/clebma/faclog(500)                                         
      cleb=0.0                                                          
      if(id+ie-if) 7000,105,7000                                        
  105 k1=ia+ib+ic                                                       
      if((-1)**k1) 7000,110,110                                         
  110 k1=ia+ib-ic                                                       
      k2=ic-iabs(ia-ib)                                                 
      k3=min0(k1,k2)                                                    
      if(k3) 7000,130,130                                               
  130 if((-1)**(ib+ie)) 7000,7000,140                                   
  140 if((-1)**(ic+if)) 7000,7000,150                                   
  150 if(ia-iabs (id)) 7000,152,152                                     
  152 if(ib-iabs (ie)) 7000,154,154                                     
  154 if(ic-iabs (if)) 7000,160,160                                     
  160 if(ia) 7000,175,165                                               
  165 if(ib) 7000,175,170                                               
  170 if(ic) 7000,180,250                                               
  175 cleb=1.0                                                          
      go to 7000                                                        
  180 fb=float(ib+1)                                                    
      cleb=((-1.0)**((ia-id)/2))/sqrt(fb)                               
      go to 7000                                                        
  250 fc2=ic+1                                                          
      iabcp=(ia+ib+ic)/2+1                                              
      iabc=iabcp-ic                                                     
      icab=iabcp-ib                                                     
      ibca=iabcp-ia                                                     
      iapd=(ia+id)/2+1                                                  
      iamd=iapd-id                                                      
      ibpe=(ib+ie)/2+1                                                  
      ibme=ibpe-ie                                                      
      icpf=(ic+if)/2+1                                                  
      icmf=icpf-if                                                      
      sqfclg=0.5*(log(fc2)-faclog(iabcp+1)                              
     1      +faclog(iabc)+faclog(icab)+faclog(ibca)                     
     2      +faclog(iapd)+faclog(iamd)+faclog(ibpe)                     
     3      +faclog(ibme)+faclog(icpf)+faclog(icmf))                    
      nzmic2=(ib-ic-id)/2                                               
      nzmic3=(ia-ic+ie)/2                                               
      nzmi= max0 (0,nzmic2,nzmic3)+1                                    
      nzmx= min0 (iabc,iamd,ibpe)                                       
      if(nzmx.lt.nzmi) go to 7000                                       
      s1=(-1.0)**(nzmi-1)                                               
      do 400 nz=nzmi,nzmx                                               
      nzm1=nz-1                                                         
      nzt1=iabc-nzm1                                                    
      nzt2=iamd-nzm1                                                    
      nzt3=ibpe-nzm1                                                    
      nzt4=nz-nzmic2                                                    
      nzt5=nz-nzmic3                                                    
      termlg=sqfclg-faclog(nz)-faclog(nzt1)-faclog(nzt2)                
     1           -faclog(nzt3)-faclog(nzt4)-faclog(nzt5)                
      ssterm=s1*exp (termlg)                                            
      cleb=cleb+ssterm                                                  
  400 s1=-s1                                                            
 7000 return                                                            
      end                                                               
c     ******************************************************************
      function rac(ia,ib,ic,id,ie,if)                                   
      implicit real*8(a-h,o-z)                                          
      common/clebma/faclog(500)                                         
      dimension  lt(6)                                                  
      rac=0.0                                                           
      k1=ia+ib-ie                                                       
      k2=ie-iabs (ia-ib)                                                
      k3=ic+id-ie                                                       
      k4=ie-iabs (ic-id)                                                
      k5=ia+ic-if                                                       
      k6=if-iabs (ia-ic)                                                
      k7=ib+id-if                                                       
      k8=if-iabs(ib-id)                                                 
      k9= min0 (k1,k2,k3,k4,k5,k6,k7,k8)                                
      if(k9) 7000,20,20                                                 
   20 k2=k1-2*(k1/2)                                                    
      k4=k3-2*(k3/2)                                                    
      k6=k5-2*(k5/2)                                                    
      k8=k7-2*(k7/2)                                                    
      if(max0(k2,k4,k6,k8)) 7000,25,7000                                
   25 ltmin=min0(ia,ib,ic,id,ie,if)                                     
      if(ltmin) 7000,30,150                                             
   30 lt(1)=ia                                                          
      lt(2)=ib                                                          
      lt(3)=ic                                                          
      lt(4)=id                                                          
      lt(5)=ie                                                          
      lt(6)=if                                                          
      ltmin=lt(1)                                                       
      kmin=1                                                            
      do 40 n=2,6                                                       
      if(lt(n)-ltmin) 35,40,40                                          
   35 ltmin=lt(n)                                                       
      kmin=n                                                            
   40 continue                                                          
      s1=1.0                                                            
      f1=ie                                                             
      f2=if                                                             
      go to (55,55,55,55,45,50),kmin                                    
   45 f1=ia                                                             
      f2=ic                                                             
      s1=(-1.0)**(k5/2)                                                 
      go to 55                                                          
   50 f1=ia                                                             
      f2=ib                                                             
      s1=(-1.0)**(k1/2)                                                 
   55 rac=s1/sqrt((f1+1.0)*(f2+1.0))                                    
      go to 7000                                                        
  150 iabep=(ia+ib+ie)/2+1                                              
      icdep=(ic+id+ie)/2+1                                              
      iacfp=(ia+ic+if)/2+1                                              
      ibdfp=(ib+id+if)/2+1                                              
      iabe=iabep-ie                                                     
      ieab=iabep-ib                                                     
      ibea=iabep-ia                                                     
      icde=icdep-ie                                                     
      iecd=icdep-id                                                     
      idec=icdep-ic                                                     
      iacf=iacfp-if                                                     
      ifac=iacfp-ic                                                     
      icfa=iacfp-ia                                                     
      ibdf=ibdfp-if                                                     
      ifbd=ibdfp-id                                                     
      idfb=ibdfp-ib                                                     
      iabcd1=(ia+ib+ic+id+4)/2                                          
      iefmad=(ie+if-ia-id)/2                                            
      iefmbc=(ie+if-ib-ic)/2                                            
      nzmax=min0(iabe,icde,iacf,ibdf)                                   
      nzmi1=-iefmad                                                     
      nzmi2=-iefmbc                                                     
      nzmin=max0(0,nzmi1,nzmi2)+1                                       
      if(nzmax.lt.nzmin) go to 7000                                     
      sqlog=faclog(iabe)+faclog(ieab)+faclog(ibea)+faclog(icde)+faclog(i
     1ecd)+faclog(idec)+faclog(iacf)+faclog(ifac)+faclog(icfa)+faclog(ib
     2df)+faclog(ifbd)+faclog(idfb)-faclog(iabep+1)-faclog(icdep+1)-facl
     3og(iacfp+1)-faclog(ibdfp+1)                                       
      sqlog=0.5*sqlog                                                   
      do 200 nz=nzmin,nzmax                                             
      nzm1=nz-1                                                         
      k1=iabcd1-nzm1                                                    
      k2=iabe-nzm1                                                      
      k3=icde-nzm1                                                      
      k4=iacf-nzm1                                                      
      k5=ibdf-nzm1                                                      
      k6=nz                                                             
      k7=iefmad+nz                                                      
      k8=iefmbc+nz                                                      
      sslog=sqlog+faclog(k1)-faclog(k2)-faclog(k3)-faclog(k4)           
     1           -faclog(k5)-faclog(k6)-faclog(k7)-faclog(k8)           
      ssterm=((-1.0)**nzm1)*exp(sslog)                                  
      rac=rac+ssterm                                                    
  200 continue                                                          
 7000 return                                                            
      end                                                               
c     ******************************************************************
      function  u9(ja,jb,jc,jd,je,jf,jg,jh,ji)                          
      implicit real*8(a-h,o-z)                                          
      u9=0.0                                                            
      k1=iabs(jb-jg)                                                    
      k2=iabs(jc-je)                                                    
      k3=iabs(jd-ji)                                                    
      minrda=max0(k1,k2,k3)                                             
      k1=jb+jg                                                          
      k2=jc+je                                                          
      k3=jd+ji                                                          
      maxrda=min0(k1,k2,k3)                                             
      if(minrda-maxrda) 30,30,20                                        
   20 go to 7000                                                        
   30 do 50 n1=minrda,maxrda,2                                          
      ramda2=n1                                                         
      u9=u9+(ramda2+1.0)*rac(jb,je,jg,jc,ja,n1)*rac(jb,jd,jg,ji,jh,n1)  
     1                  *rac(jc,jd,je,ji,jf,n1)                         
   50 continue                                                          
 7000 return                                                            
      end                                                               
c     ******************************************************************
      subroutine coulfn (eta,rhoz,drho,lmax)                            
      implicit real*8 ( a-h , o-z )                                     
      common/jat/ sa(90),fa(90),fpp(90),ga(90),gpp(90)                  
      dimension gd(5)                                                   
      data cg0, cg1, cg2, cg3, cg4, cg5/                                
     +1.223404016d0,4.959570165d-2,8.888888889d-3,                      
     +2.455199181d-3,9.108958061d-4,2.534684115d-4/                     
      data cgp,cp1,cp2,cp3,cp4,cp5 /                                    
     +-7.078817734d-1,1.728260369d-1,3.174603174d-4,                    
     +3.581214850d-3,3.117824680d-4,9.073966427d-4/                     
      data gd(1),gd(2),gd(3),gd(4),gd(5) / 12.,-360.,1260.,-1680.,1188./
      nr=1                                                              
      epsi=0.d0                                                         
      eps=0.1*epsi                                                      
      eps=max (eps, 1.0d-11)                                            
      eta2=eta+eta                                                      
      etas=eta*eta                                                      
      lp=max0 (lmax,12)+1                                               
      t=lp                                                              
      u=t*t+etas                                                        
      v=t/u                                                             
      w=eta/u                                                           
      x=v*v-w*w                                                         
      y=2.*v*w                                                          
      u=sqrt (u)                                                        
      sig0=eta*(dlog(u)-1.)+(t-0.5)* datan(eta/t)                       
      do 20 i=1,5                                                       
      sig0=sig0-w/gd(i)                                                 
      t=v*x-w*y                                                         
      w=v*y+w*x                                                         
      v=t                                                               
   20 continue                                                          
   30 if (lp .le. lmax+1) sa(lp)=sig0                                   
      lp=lp-1                                                           
      if (lp .le. 0) go to 100                                          
      t=lp                                                              
      sig0=sig0- datan (eta/t)                                          
      go to 30                                                          
  100 emax=(1.0e-5/eps)**0.16666667                                     
      if (eta .lt. emax) go to 200                                      
      r=eta2                                                            
      t=6.0                                                             
      t=eta**(1.0/t)                                                    
      w=eta*t                                                           
      u=t-t*(cg2+cg4/etas)/etas                                         
      v=(cg1+(cg3+cg5/etas)/etas)/w                                     
      g=cg0*(u+v)                                                       
      t=1./t                                                            
      w=eta*t                                                           
      u=t+t*(cp2+cp4/etas)/etas                                         
      v=(cp1+(cp3+cp5/etas)/etas)/w                                     
      gp=cgp*(u-v)                                                      
      go to 300                                                         
  200 r=max0 (nr,1)-1                                                   
      r=rhoz+r*drho                                                     
      t=12.+1.4*eta                                                     
      if (t .gt. r) r=t                                                 
      fk=1.                                                             
      f=1.                                                              
      gk=0.                                                             
      g=0.                                                              
      fsk=0.                                                            
      fp=0.                                                             
      gsk=1.-eta/r                                                      
      gp=gsk                                                            
      epss=eps*eps                                                      
      n=r+r                                                             
      do 210 kp=1,n                                                     
      t=kp+kp                                                           
      u=t*r                                                             
      ak=(t-1.)*eta/u                                                   
      v=kp*(kp-1)                                                       
      bk=(etas-v)/u                                                     
      t=ak*fk-bk*gk                                                     
      gk=ak*gk+bk*fk                                                    
      fk=t                                                              
      t=ak*fsk-bk*gsk-fk/r                                              
      gsk=ak*gsk+bk*fsk-gk/r                                            
      fsk=t                                                             
      f=f+fk                                                            
      g=g+gk                                                            
      fp=fp+fsk                                                         
      gp=gp+gsk                                                         
      test=fk*fk+gk*gk+fsk*fsk + gsk*gsk                                
      if (test .lt. epss) go to 220                                     
  210 continue                                                          
  220 t=r-eta*dlog(r+r)+sig0                                            
      u=cos (t)                                                         
      v=sin (t)                                                         
      g=f*u-g*v                                                         
      gp=fp*u-gp*v                                                      
  300 rs=r                                                              
      rho=rhoz                                                          
      f=g                                                               
      fp=gp                                                             
      is=0                                                              
      ir=1                                                              
      t=r-rho                                                           
      if (t) 600,700,310                                                
  310 if (nr .le. 1) go to 320                                          
      is=t/drho                                                         
      is=min0 (is+1,nr)                                                 
  320 t=is                                                              
      rho=rhoz+t*drho                                                   
      gg=rho                                                            
      is=is+1                                                           
      ir=is                                                             
  330 rho=rho-drho                                                      
      ir=ir-1                                                           
      if (ir .gt. 0) go to 600                                          
      ir=max0 (is,1)                                                    
      r=rs                                                              
      rho=gg                                                            
      f=g                                                               
      fp=gp                                                             
      go to 350                                                         
  340 rho=rho+drho                                                      
      ir=ir+1                                                           
  350 if (ir .gt. nr) return                                            
  600 h=0.5                                                             
      w=r-eta2                                                          
      if (r-1.0) 601,602,602                                            
  601 h=0.5*r                                                           
  602 if (w) 603,605,605                                                
  603 t=sqrt (-r/(w+w))                                                 
      if (t-h) 604,605,605                                              
  604 h=t                                                               
  605 last=0                                                            
      t=rho-r                                                           
      if (t) 606,700,607                                                
  606 h=-h                                                              
  607 u=t-h                                                             
      if (u*h) 608,608,609                                              
  608 h=t                                                               
      last=1                                                            
  609 u=0.0                                                             
      t=1.0                                                             
      b1=0.0                                                            
      b2=f                                                              
      b3=h*fp                                                           
      f=f+b3                                                            
      v=0.0                                                             
  610 it=0                                                              
  620 v=-h*(h*b1+w*b2+u*v)/(r*t)                                        
      fp=fp+v                                                           
      u=t                                                               
      t=t+1.0                                                           
      b1=b2                                                             
      b2=b3                                                             
      b3=h*v/t                                                          
      f=f+b3                                                            
      test=b3                                                           
      testp=v                                                           
      if (w) 630,640,640                                                
  630 test=b3/f                                                         
      testp=v/fp                                                        
  640 if (abs(test)+abs(testp)-eps) 650,610,610                         
  650 if (it) 660,660,670                                               
  660 it=1                                                              
      go to 620                                                         
  670 r=r+h                                                             
      if (last) 600,600,700                                             
  700 k=lmax+1                                                          
      x=f                                                               
      y=fp                                                              
      do 710 j=1,k                                                      
      ga(j)=x                                                           
      al=j                                                              
      t=j*j                                                             
      u=t/rho+eta                                                       
      v=sqrt (t+etas)                                                   
      w=(u*x-al*y)/v                                                    
      y=(v*x-u*w)/al                                                    
      x=w                                                               
  710 continue                                                          
      lp=rho                                                            
      lp=max0 (lp+10, lmax+20)                                          
      b3=0.                                                             
      b2=1.0e-20                                                        
      w=1.0/rho                                                         
      al=lp+1                                                           
      v=eta/al                                                          
      u=0.                                                              
      do 840 j=1,lp                                                     
      k=lp+1-j                                                          
      al=k                                                              
      t=k+k+1                                                           
      b1=t*(v/al+w)*b2-u*b3                                             
      v=eta/al                                                          
      u=sqrt (1.0+v*v)                                                  
      b1=b1/u                                                           
      b3=b2                                                             
      b2=b1                                                             
      if (k-lmax-1) 810,810,820                                         
  810 fa(k)=b1                                                          
      go to 840                                                         
  820 test=b1                                                           
      if (abs(test)-1.) 840,840,830                                     
  830 b2=b2*1.0d-20                                                     
      b3=b3*1.0d-20                                                     
  840 continue                                                          
      t=(w+eta)*b2-u*b3                                                 
      u=1./(t*f-b1*fp)                                                  
      k=lmax+1                                                          
      do 850 j=1,k                                                      
      fa(j)=u*fa(j)                                                     
  850 continue                                                          
      do 400 l=1,k                                                      
      fl=l                                                              
      flsq=fl*fl                                                        
      fac1=eta/fl+fl/rhoz                                               
      fac2=dsqrt(eta*eta+flsq)/fl                                       
      fpp(l)=fac1*fa(l)-fac2*fa(l+1)                                    
      gpp(l)=fac1*ga(l)-fac2*ga(l+1)                                    
  400 continue                                                          
      if (ir-is) 330,340,340                                            
      end                                                               
c     ******************************************************************
      subroutine ffgenz(i,k,ich,istep)                                  
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/ffc/amp,d02,fts,ta,tb,fnrng,ktff(5),nrange,nafac,nbfac,lsa 
     1,lsb,iread(15),jjj,fread(28),wread(21),ktrl(20),ff(905),ffi(905), 
     2ffone(905),nrangd(2,4),n1xd(2,4),itap(66,10,2),frfac(2,4),frso(2) 
      complex*16 ffone,frfac                                            
      common/zrff/fis,fls,fjs,ffout(10),bb,bsubl,nrn,lss,redwd,zeta,    
     1fma,fmas,fmb,fmbs,sa,fja,sb,fjb,za,zas,zb,zbs,vb(905)             
      fja=tspn(i)                                                       
      sa=pspn(i)                                                        
      fma=tmas(i)                                                       
      fmas=pmas(i)                                                      
      za=tz(i)                                                          
      zas=pz(i)                                                         
      fjb=tspn(k)                                                       
      sb=pspn(k)                                                        
      fmb=tmas(k)                                                       
      fmbs=pmas(k)                                                      
      zb=tz(k)                                                          
      zbs=pz(k)                                                         
      fis=trs(istep,ich)                                                
      fls=ltr(istep,ich)                                                
      fjs=trj(istep,ich)                                                
      lss=fls+0.01                                                      
      if(abs(fmas-fmbs).lt.0.01) go to 30                               
      if(fma.gt.fmb) go to 20                                           
      fac=fmb*tmas(1)/(fma*tmas(i))                                     
      sfac=sqrt((2.0*sa+1.0)/(2.0*fis+1.0))                             
      phase=1.0                                                         
      dr=drd(i)                                                         
      go to 40                                                          
   20 fac=fma*tmas(1)/(fmb*tmas(k))                                     
      sfac=sqrt((2.0*fja+1.0)*(2.0*sb+1.0)/((2.0*fjb+1.0)*(2.0*fis+1.0))
     1)                                                                 
      m=abs(fja+fjs-fjb+sb+fis-sa)+0.01                                 
      phase=(-1.0)**m                                                   
      dr=drd(k)                                                         
      go to 40                                                          
   30 fac=1.0                                                           
      sfac=sqrt(2.0*sa+1.0)                                             
      phase=1.0                                                         
      dr=drd(ich)                                                       
   40 fac=fac*sfac*phase                                                
      do 45 nr=1,nrmax                                                  
      ff(nr)=0.0                                                        
   45 ffi(nr)=0.0                                                       
      call inputb                                                       
      ktf=ktff(1)                                                       
      if(ktf.eq.5) go to 767                                            
      call ffsub4                                                       
*     ---------------------------------------------------------
*     signs of radial wave functions meas the origin
*     ---------------------------------------------------------
      print*,' wave function sign near origin ',ff(3)/abs(ff(3))
*     print*,ff(3),ffi(3)
      print*
  767 continue
*     write(6,100) ich,istep                                            
  100 format(1h ,/,5x,30h*****(   form factor  channel=,i2,2x           
     1,5hstep=,i2,3x,6h)*****)                                          
      call writeb                                                       
      call fnloc(i,k)                                                   
      fac1=sqrt(d02)*amp                                                
      if(abs(fmas-fmbs).lt.0.01) fac1=amp                               
      if(ktf.eq.8) fac1=1.0                                             
      fac=fac*fac1                                                      
      do 60 nr=1,nrmax                                                  
      ff(nr)=ff(nr)*fac                                                 
   60 ffi(nr)=ffi(nr)*fac                                               
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine fnloc(i,k)                                             
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/ffc/amp,d02,fts,ta,tb,fnrng,ktff(5),nrange,nafac,nbfac,lsa 
     1,lsb,iread(15),jjj,fread(28),wread(21),ktrl(20),ff(905),ffi(905), 
     2ffone(905),nrangd(2,4),n1xd(2,4),itap(66,10,2),frfac(2,4),frso(2) 
      complex*16 ffone,frfac                                            
      common/zrff/fis,fls,fjs,ffout(10),bb,bsubl,nrn,lss,redwd,zeta,    
     1fma,fmas,fmb,fmbs,sa,fja,sb,fjb,za,zas,zb,zbs,vb(905)             
      common/pot/vd(2),wd(2),vsod(2),wsod(2),rrd(2),ard(2),rcd(2),      
     1csdgd(2),rid(2),aid(2),rsord(2),asord(2),rsoid(2),asoid(2),       
     2ktop(2),rgd(2),agd(2),ktisp(2),risrd(2),aisrd(2),visd(2),risid(2),
     3aisid(2),wisd(2),csd,drhod(12),rhobn,rhoin,rhosor,rhosoi,rhoisr,  
     4rhoisi,rhobc,rhobng,fkayar,fkayai,fkasor,fkasoi,fkaisr,fkaisi,    
     5fkaybg,vbe,wbe,vsoe,wsoe,vise,wise,vbcd(905,2),wbcd(905,2),       
     6vsov(905,2),wsov(905,2),pnloc(2)                                  
*     data facm/0.048196758/       
      data facm/0.0478450d0/                                          
      if(abs(fnrng).lt.1.e-20) goto 398                                 
      if(fmas-fmbs)143,398,142                                          
  142 i1=i                                                              
      i2=k                                                              
      go to 144                                                         
  143 i1=k                                                              
      i2=i                                                              
  144 continue                                                          
      em1=ecmd(i1)                                                      
      em2=ecmd(i2)                                                      
      temp=pmas(i2)*fnrng**2/pmas(i1)                                   
      fact=temp*tmas(i2)/tmas(i1)                                       
      fmx=pmas(i1)-pmas(i2)                                             
      temp=temp*fmx*facm                                                
      go to 350                                                         
  340 continue                                                          
      do 345 m=1,nrmax                                                  
      ctemp1=temp*(vbcd(m,i2)*em2-vbcd(m,i1)*em1)                       
      ctemp2=temp*(wbcd(m,i2)*em2-wbcd(m,i1)*em1)                       
      fact=exp(-ctemp1)                                                 
      ut1= fact*cos(ctemp2)                                             
      ut2=-fact*sin(ctemp2)                                             
      uf1=ut1*ff (m)-ut2*ffi(m)                                         
      uf2   =ut1*ffi(m)+ut2*ff (m)                                      
      ff (m)=uf1                                                        
      ffi(m)=uf2                                                        
  345 continue                                                          
      go to 398                                                         
  350 continue                                                          
      if(fnrng.gt.0.0) go to 360                                        
      do 355 m=1,nrmax                                                  
      fact=exp(-temp*vb(m))                                             
      ff (m)=ff (m)*fact                                                
      ffi(m)=ffi(m)*fact                                                
  355 continue                                                          
      go to 340                                                         
  360 continue                                                          
      do 365 m=1,nrmax                                                  
      ctemp1=vb(m)+vbcd(m,i2)*em2-vbcd(m,i1)*em1                        
      ctemp2=wbcd(m,i2)*em2-wbcd(m,i1)*em1                              
c     ctemp1=1.0-temp*ctemp1                                            
      ctemp1=1.0+temp*ctemp1                                            
c     ctemp2=-temp*ctemp2                                               
      ctemp2=temp*ctemp2                                                
c     det=1.d0                                                          
      det=ctemp1**2+ctemp2**2                                           
c     ut1=(ctemp1*ff(m)-ctemp2*ffi(m))/det                              
      ut1=(ctemp1*ff(m)+ctemp2*ffi(m))/det                              
c     ut2=(ctemp1*ffi(m)+ctemp2*ff(m))/det                              
      ut2=(ctemp1*ffi(m)-ctemp2*ff(m))/det                              
      ff(m)=ut1                                                         
      ffi(m)=ut2                                                        
  365 continue                                                          
  398 continue                                                          
      fact1=facm*fmud(i)*ecmd(i)*pnloc(i)**2/8.0                        
      fact2=facm*fmud(k)*ecmd(k)*pnloc(k)**2/8.0                        
      if(abs(fact1+fact2).lt.1.e-20) goto 410                           
      do 400 m=1,nrmax                                                  
      ctemp1=-fact1*(vbcd(m,i)-1.0)-fact2*(vbcd(m,k)-1.0)               
      ctemp2=-fact1*wbcd(m,i)-fact2*wbcd(m,k)                           
      fact=exp(ctemp1)                                                  
      ut1=fact*cos(ctemp2)                                              
      ut2=fact*sin(ctemp2)                                              
      uf1=ut1*ff(m)-ut2*ffi(m)                                          
      uf2=ut1*ffi(m)+ut2*ff(m)                                          
      ff (m)=uf1                                                        
      ffi(m)=uf2                                                        
  400 continue                                                          
  410 continue                                                          
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine ffsub4                                                 
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drm(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/ffc/amp,d02,fts,ta,tb,fnrng,ktff(5),nrange,nafac,nbfac,lsa 
     1,lsb,iread(15),jjj,fread(28),wread(21),ktrl(20),ff(905),ffi(905), 
     2ffone(905),nrangd(2,4),n1xd(2,4),itap(66,10,2),frfac(2,4),frso(2) 
      complex*16 ffone,frfac                                            
      common/zrff/fis,fls,fjs,ffout(10),bb,bsubl,nrn,lss,redwd,zeta,    
     1fma,fmas,fmb,fmbs,sa,fja,sb,fjb,za,zas,zb,zbs,vb(905)             
      dimension drd(1),wfc(905),ffr(905,2),ast(20),wfcr(905),wfsr(905)  
      dimension rar2(905),rarr(905)
      ktrlf4=0                                                          
      jhw=1                                                             
      eps7=0.0001                                                       
      fnrmax=nrmax                                                      
      drd(1)=dr                                                         
      lmom=fls+0.01                                                     
      nramax=nrmax                                                      
      ramax=dr*fnrmax                                                   
      nod=iread(1)                                                      
      wzz=iread(2)                                                      
      bengy=fread(1)                                                    
      wsmas=fread(2)                                                    
      wmass=fread(3)                                                    
      dmat1=fread(4)                                                    
      drz1=fread(5)                                                     
      wr0=fread(6)                                                      
      wrc=fread(7)                                                      
      wa=fread(8)                                                       
      wls=fread(9)*2.0                                                  
      pnlocb=fread(23)  
      if(frso(1).lt.0.d0.and.frso(2).lt.0.d0) then
       frso(1)=wr0
       frso(2)=wa 
      endif  
      wrso=frso(1)
      waso=frso(2)                                                
      lmom1=lmom+1                                                      
      test=1.0e16                                                       
      lmom2=lmom1+1                                                     
      fnod=nod                                                          
      flmom=lmom                                                        
      flmom1=lmom1                                                      
      t=wmass**0.333333333                                              
      radi=wr0*t  
      radso=wrso*t                                                      
      radc=wrc*t                                                        
      radz=radi                                                         
      wr=drd(jhw)                                                       
      nr1=nramax+1                                                      
      do 8 i=1,nr1                                                      
      ex=exp((wr-radi)/wa)   
      exso=exp((wr-radso)/waso)                                         
      wfcr(i)=1.0/(1.0+ex)  
      wfso=1.0/(1.0+exso)                                            
      wfsr(i)=exso*wfso*wfso                                        
      if(wr-radc) 3,3,4                                                 
*   3 wfc(i)=0.7199262*wzz*(3.0-wr*wr/(radc*radc))/radc    
    3 wfc(i)=1.43997d0/2.d0*wzz*(3.0-wr*wr/(radc*radc))/radc            
      go to 5                                                           
*   4 wfc(i)=1.4398523*wzz/wr          
    4 wfc(i)=1.43997d0*wzz/wr                                          
    5 wr=wr+drd(jhw)                                                    
    8 continue                                                          
      do 2000 ii=1,15                                                   
      do 1000 jj=1,15                                                   
      drz=drz1+0.3*float(ii/2)*(-1.0)**ii                               
      dmat=dmat1+2.01*float(jj/2)*(-1.0)**jj                            
      korec=0                                                           
      niter=0                                                           
      incr=0                                                            
      if(ktrlf4-1) 10,20,10                                             
   10 vdepth=bengy+(3.1415926*(fnod+0.5*flmom1))**2/(0.048228*wsmas*(rad
     1i+drz)**2)                                                        
      go to 30                                                          
   20 bengy=vdepth-(3.1415926*(fnod+0.5*flmom1))**2/(0.048228*wsmas     
     1*radz*radz)                                                       
      if(bengy-eps7) 25,25,30                                           
   25 if(abs(drz).lt.0.01) drz=0.3                                      
      radz=radz+drz                                                     
      incr=incr+1                                                       
      if(incr-20) 20,20,27                                              
   27 kcheck=11                                                         
      go to 7400                                                        
   30 flns=fjs*(fjs+1.0)-fls*(fls+1.0)-fis*(fis+1.0)                    
      match=radi/drd(jhw)+dmat                                          
   80 fmut=wsmas*wmass/(wsmas+wmass)                                    
*     wk=0.2195376*sqrt(fmut*bengy)     
      wk=0.218735d0*sqrt(fmut*bengy)                                
      wrhoc=wk*radc                                                     
      wrhoz=wk*radz                                                     
      wrhocs=wrhoc*wrhoc                                                
*     weta=0.7199262*wzz*wk/bengy   
      weta=1.43997d0/2.d0*wzz*wk/bengy                                  
      wetac=weta/wrhoc                                                  
      wdrho=drd(jhw)*wk                                                 
      wvs=wls*wk/(bengy*waso)                                           
      drhosq=wdrho*wdrho                                                
      dr56=0.8333333333*drhosq                                          
      dr12=0.1*dr56                                                     
      fl1=lmom*lmom1                                                    
  100 wrho=wdrho                                                        
      wvc=vdepth/bengy                                                  
      zer=1.0                                                           
      do 180 j=1,lmom2                                                  
      a1=-wvs*flns*wfsr(j)/(flmom1+flmom1)                              
      b1=1.0-wvc*wfcr(j)+3.0*wetac                                      
      b2=wvs*flns*wfsr(j)                                               
      a2=(b1-b2*a1)/(4.0*flmom1+2.0)                                    
      a3=(b1*a1-b2*a2)/(6.0*flmom1+6.0)                                 
      wrhosq=wrho*wrho                                                  
      b3=weta/(wrhoc*wrhocs)                                            
      a4=(b1*a2-b2*a3-b3)/(8.0*flmom1+12.0)                             
      a5=(b1*a3-b2*a4-b3*a1)/(10.0*flmom1+20.0)                         
      a6=(b1*a4-b2*a5-b3*a2)/(12.0*flmom1+30.0)                         
      ffr(j,jhw)=(wrho**lmom1)*(1.0+a1*wrho+a2*wrho*wrho+a3*wrhosq*wrho 
     1+a4*wrhosq*wrhosq+a5*wrho*wrhosq*wrhosq+a6*wrhosq**3)             
  180 wrho=wrho+wdrho                                                   
      mat1=match+1                                                      
      x1=wdrho*flmom1                                                   
      x2=x1+wdrho                                                       
      x3=x2+wdrho                                                       
      do 200 i=lmom1,mat1                                               
      fac1=1.0-dr12*(fl1/(x1*x1)+1.0-wvc*wfcr(i)-wvs*flns*wfsr(i)/x1    
     1+wfc(i)/bengy)                                                    
      fac2=2.0+dr56*(fl1/(x2*x2)+1.0-wvc*wfcr(i+1)-wvs*flns*wfsr(i+1)/x2
     1+wfc(i+1)/bengy)                                                  
      fac3=1.0-dr12*(fl1/(x3*x3)+1.0-wvc*wfcr(i+2)-wvs*flns*wfsr(i+2)/x3
     1+wfc(i+2)/bengy)                                                  
      ffr(i+2,jhw)=(ffr(i+1,jhw)*fac2-ffr(i,jhw)*fac1)/fac3             
      if(abs(ffr(i+2,jhw))-test) 195,185,185                            
  185 ip2=i+2                                                           
      do 190 ip=1,ip2                                                   
      ffr(ip,jhw)=ffr(ip,jhw)/test                                      
  190 continue                                                          
      zer=zer/test                                                      
  195 x1=x2                                                             
      x2=x3                                                             
      x3=x3+wdrho                                                       
  200 continue                                                          
      nodes=0                                                           
      do 220 i=2,match                                                  
      if(ffr(i,jhw)*ffr(i+1,jhw)) 210,215,220                           
  210 nodes=nodes+2                                                     
      go to 220                                                         
  215 nodes=nodes+1                                                     
  220 continue                                                          
      nn=nodes/2                                                        
      fnn=nn                                                            
      if(nod-nn) 225,240,225                                            
  225 korec=korec+1                                                     
      if(korec-8) 228,228,226                                           
  226 kcheck=10                                                         
      go to 2000                                                        
  228 vcor=(wrhoz*wrhoz+9.86959*(fnod+0.5*flmom1)**2)/(wrhoz*wrhoz      
     1+9.86959*(fnn+0.5*flmom1)**2)                                     
      if(ktrlf4-1) 230,235,235                                          
  230 vdepth=vcor*vdepth                                                
      go to 100                                                         
  235 bengy=bengy/vcor                                                  
      go to 80                                                          
  240 dffr1=((ffr(match+3,jhw)-ffr(match-3,jhw))/60.0+3.0*(ffr(match-2, 
     1jhw)-ffr(match+2,jhw))/20.0+3.0*(ffr(match+1,jhw)-ffr(match-1,jhw)
     2)/4.0)/wdrho                                                      
      ast(1)=1.0                                                        
      t1=1.0                                                            
      t2=2.0                                                            
      do 320 i=1,15                                                     
      ast(i+1)=(flmom1-t1-weta)*(flmom1+t1-1.0+weta)*ast(i)/t2          
      t1=t1+1.0                                                         
      t2=t2+2.0                                                         
  320 continue                                                          
      rhoa=ramax*wk                                                     
      wrho=rhoa                                                         
      jrho=nramax                                                       
  325 frt=1.0                                                           
      trhoa=wrho                                                        
      do 330 i=1,15                                                     
      frt=frt+ast(i+1)/trhoa                                            
      trhoa=trhoa*wrho                                                  
  330 continue                                                          
      ffr(jrho,jhw+1)=frt/exp(wrho+weta*log(2.0*wrho))                  
      if(jrho-nramax) 340,340,350                                       
  340 jrho=jrho+1                                                       
      wrho=wrho+wdrho                                                   
      go to 325                                                         
  350 x1=rhoa-wdrho                                                     
      x2=rhoa                                                           
      x3=x2+wdrho                                                       
      imax=nramax-match+3                                               
      do 360 i=1,imax                                                   
      k=nramax-i                                                        
      fac1=1.0-dr12*(fl1/(x1*x1)+1.0-wvc*wfcr(k)-wvs*flns*wfsr(k)/x1    
     1+wfc(k)/bengy)                                                    
      fac2=2.0+dr56*(fl1/(x2*x2)+1.0-wvc*wfcr(k+1)-wvs*flns*wfsr(k+1)/x2
     1+wfc(k+1)/bengy)                                                  
      fac3=1.0-dr12*(fl1/(x3*x3)+1.0-wvc*wfcr(k+2)-wvs*flns*wfsr(k+2)/x3
     1+wfc(k+2)/bengy)                                                  
      ffr(k,jhw+1)=(ffr(k+1,jhw+1)*fac2-ffr(k+2,jhw+1)*fac3)/fac1       
      if(abs(ffr(k,jhw+1))-test) 358,352,352                            
  352 nrten=nramax-k+2                                                  
      do 356 iten=1,nrten                                               
      kten=iten+k-1                                                     
      ffr(kten,jhw+1)=ffr(kten,jhw+1)/test                              
  356 continue                                                          
  358  x3=x2                                                            
      x2=x1                                                             
      x1=x1-wdrho                                                       
  360 continue                                                          
      dffr2=((ffr(match+3,jhw+1)-ffr(match-3,jhw+1))/60.0+3.0*(ffr(match
     1-2,jhw+1)-ffr(match+2,jhw+1))/20.0+3.0*(ffr(match+1,jhw+1)        
     2-ffr(match-1,jhw+1))/4.0)/wdrho                                   
      ratio=ffr(match,jhw)/ffr(match,jhw+1)                             
      tlogd1=dffr1/ffr(match,jhw)                                       
      tlogd2=dffr2/ffr(match,jhw+1)                                     
      difnce=abs(tlogd1-tlogd2)                                         
      if(difnce-eps7) 510,510,400                                       
  400 niter=niter+1                                                     
      if(niter-20) 410,410,405                                          
  405 kcheck=12                                                         
      go to 1000                                                        
  410 fnum=ffr(match,jhw+1)*dffr2*ratio*ratio-ffr(match,jhw)*dffr1      
      sum=0.0                                                           
      do 480 i=1,nramax,2                                               
      if(i-1) 420,420,430                                               
  420 sum1=0.0                                                          
      go to 445                                                         
  430 if(i-match) 440,440,450                                           
  440 sum1=ffr(i-1,jhw)*ffr(i-1,jhw)                                    
  445 sum2=ffr(i,jhw)*ffr(i,jhw)                                        
      sum3=ffr(i+1,jhw)*ffr(i+1,jhw)                                    
      go to 455                                                         
  450 sum1=ffr(i-1,jhw+1)*ffr(i-1,jhw+1)*ratio*ratio                    
      sum2=ffr(i,jhw+1)*ffr(i,jhw+1)*ratio*ratio                        
      sum3=ffr(i+1,jhw+1)*ffr(i+1,jhw+1)*ratio*ratio                    
  455 if(ktrlf4-1) 460,470,470                                          
  460 if(i-1) 462,462,465                                               
  462 sum1=0.0                                                          
      go to 467                                                         
  465 sum1=-sum1*wfcr(i-1)*wvc                                          
  467 sum2=-sum2*wfcr(i)*wvc                                            
      sum3=-sum3*wfcr(i+1)*wvc                                          
  470 sum=sum+sum1+4.0*sum2+sum3                                        
  480 continue                                                          
      denom=sum*wdrho/3.0                                               
      incr=0                                                            
      ram1=fnum/denom                                                   
  482 ramda=1.0+ram1                                                    
      if(ramda-eps7) 485,485,488                                        
  485 ram1=0.5*ram1                                                     
      incr=incr+1                                                       
      if(incr-10) 482,482,486                                           
  486 kcheck=13                                                         
      go to 1000                                                        
  488 korec=0                                                           
      if(ktrlf4-1) 490,500,500                                          
  490 vdepth=ramda*vdepth                                               
      go to 100                                                         
  500 bengy=bengy*ramda                                                 
      go to 80                                                          
 1000 continue                                                          
 2000 continue                                                          
      go to 7400                                                        
  510 kcheck=0                                                          
      do 520 i=match,nr1                                                
      ffr(i,jhw)=ratio*ffr(i,jhw+1)                                     
  520 continue                                                          
      sum=0.0                                                           
      do 570 i=1,nramax,2                                               
      if(i-1) 540,540,550                                               
  540 sum1=0.0                                                          
      go to 560                                                         
  550 sum1=ffr(i-1,jhw)*ffr(i-1,jhw)                                    
  560 sum2=ffr(i,jhw)*ffr(i,jhw)                                        
      sum3=ffr(i+1,jhw)*ffr(i+1,jhw)                                    
      sum=sum+sum1+4.0*sum2+sum3                                        
  570 continue                                                          
      sum=sum*drd(jhw)/3.0                                              
      znorm=1.0/sqrt(sum)                                               
      r=0.0                                              
      do 600 i=1,nr1 
      r=r+drd(jhw)                                             
      ff(i)=znorm*ffr(i,jhw)/r                                          
      vb(i)=-bengy+vdepth*wfcr(i)+wls*flns*wfsr(i)/(waso*r)-wfc(i)      
  600 continue                                                          
*     fact=0.048196758*fmut*pnlocb**2/8.0  
      fact=0.0478845d0*fmut*pnlocb**2/8.0                             
      if(fact.lt.1.e-10) go to 700                                      
      sum=0.0                                                           
      r=0.0                                                             
      do 650 i=1,nr1                                                    
      r=r+drd(jhw)                                                      
      ff(i)=ff(i)*exp(-fact*(vb(i)+bengy))                              
  650 sum=sum+(ff(i)*r)**2                                              
      sum=1.0/sqrt(sum*drd(jhw))                                        
      do 660 i=1,nr1                                                    
  660 ff(i)=ff(i)*sum                                                   
  700 r=0.0 
*     calculate and print orbital rms radius       
      rarr(1)=0.d0
      rar2(1)=0.d0 
      do i=1,nrmax                                                    
       r=r+drd(jhw)      
       rarr(i+1)=r*r*ff(i)*ff(i)   
       rar2(i+1)=r*r*rarr(i+1)    
      enddo
      call sim(rarr,res0,1,nramax,drd(jhw))
      call sim(rar2,res2,1,nramax,drd(jhw))
      print*,'------------------------------------------------------'
      print*, 'Wave function normalisation is: ',real(res0)
      print*, 'Wave function rms radius is (fm)',real(sqrt(res2))
      print*,'------------------------------------------------------'
*     calculate and print orbital ANC at two radii (neutral transfer)
      r=0.d0
      rar2(1)=0.d0 
      rarr(1)=0.d0
*     following is radial function r*u(r) for ANC   
      do i=1,nrmax                                                    
       rar2(i+1)=rar2(i)+drd(jhw)  
       rarr(i+1)=rar2(i+1)*ff(i)  
      enddo
*     do i=1,nrmax+1
*      write(17,*) rar2(i),rarr(i)  
*     enddo
*     wk=0.2195376d0*sqrt(fmut*bengy)
      wk=0.218735d0*sqrt(fmut*bengy)
*     do initially at two radii for testing 
      ma1=nrmax-40
      ma2=nrmax-50
      cr1=wk*rar2(ma1)
      cr2=wk*rar2(ma2)
      if(lmom.eq.0) then
       basy1=1.d0 
       basy2=1.d0
      else if(lmom.eq.1) then
       basy1=1.d0+1.d0/cr1 
       basy2=1.d0+1.d0/cr2
      else if(lmom.eq.2) then
       basy1=1.d0+3.d0/cr1+3.d0/cr1**2 
       basy2=1.d0+3.d0/cr2+3.d0/cr2**2 
      else if(lmom.eq.3) then
       basy1=1.d0+6.d0/cr1+15.d0/cr1**2+15.d0/cr1**3 
       basy2=1.d0+6.d0/cr2+15.d0/cr2**2+15.d0/cr2**3
      else if(lmom.eq.4) then
       basy1=1.d0+10.d0/cr1+45.d0/cr1**2+105.d0/cr1**3+105.d0/cr1**4 
       basy2=1.d0+10.d0/cr2+45.d0/cr2**2+105.d0/cr2**3+105.d0/cr2**4 
      else if(lmom.eq.5) then
       basy1=1.d0+15.d0/cr1+105.d0/cr1**2+420.d0/cr1**3+
     *       945.d0/cr1**4+945.d0/cr1**5 
       basy2=1.d0+15.d0/cr2+105.d0/cr2**2+420.d0/cr2**3+
     *       945.d0/cr2**4+945.d0/cr2**5 
      endif
      basy1=rarr(ma1)/exp(-cr1)/basy1
      basy2=rarr(ma2)/exp(-cr2)/basy2
      print*, 'Wave function ANC (neutral transfers) at two radii:' 
      print*,'------------------------------------------------------'
      print*, 'Bound state Kappa = ',wk
      print*, 'Bound state step =  ',drd(jhw)
      print*,'------------------------------------------------------'
      print*, 'r=',real(rar2(ma1)),' [1] asymptotic norm**2 = ',
     *         real(basy1*basy1)
      print*, 'r=',real(rar2(ma2)),' [2] asymptotic norm**2 = ',
     *         real(basy2*basy2)
      print*,'------------------------------------------------------'
      ffout(1)=radi                                                     
      ffout(2)=float(match)*drd(jhw)                                    
      ffout(3)=ramax                                                    
      ffout(4)=korec                                                    
      ffout(5)=niter                                                    
      ffout(6)=vdepth                                                   
      ffout(7)=ramda                                                    
      ffout(8)=tlogd1                                                   
      ffout(9)=tlogd2                                                   
      ffout(10)=ratio                                                   
      go to 8000                                                        
 7400 continue
*     write(6,7500) kcheck                                              
 7500 format(//5x,7hkcheck=,i3,3x,9hin ffsub4)                          
 8000 return                                                            
      end                                                               
c     ******************************************************************
      function clebz(l1cl,l2cl,l3cl)                                    
      implicit real*8(a-h,o-z)                                          
      common/clebma/faclog(500)                                         
      lsum=(l1cl+l2cl+l3cl)/2                                           
      lg=lsum/2                                                         
      if(2*lg-lsum) 15,20,20                                            
   15 clebz=0.0                                                         
      go to 60                                                          
   20 l1=lsum-l1cl                                                      
      l2=lsum-l2cl                                                      
      l3=lsum-l3cl                                                      
      if(min0(l1,l2,l3)) 15,22,22                                       
   22 fact1=sqrt(float(l3cl+1))                                         
      lx=l3cl/2+lg                                                      
      if(2*(lx/2)-lx) 30,32,32                                          
   30 fact1=-fact1                                                      
   32 l12=l1/2+1                                                        
      l22=l2/2+1                                                        
      l32=l3/2+1                                                        
      h1=0.5*(faclog(l1+1)+faclog(l2+1)+faclog(l3+1)-faclog(lsum+2))+   
     1faclog(lg+1)-faclog(l12)-faclog(l22)-faclog(l32)                  
      clebz=fact1*exp(h1)                                               
   60 return                                                            
      end                                                               
c     ******************************************************************
      subroutine inputa                                                 
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/pot/vd(2),wd(2),vsod(2),wsod(2),rrd(2),ard(2),rcd(2),      
     1csdgd(2),rid(2),aid(2),rsord(2),asord(2),rsoid(2),asoid(2),       
     2ktop(2),rgd(2),agd(2),ktisp(2),risrd(2),aisrd(2),visd(2),risid(2),
     3aisid(2),wisd(2),csd,drhod(12),rhobn,rhoin,rhosor,rhosoi,rhoisr,  
     4rhoisi,rhobc,rhobng,fkayar,fkayai,fkasor,fkasoi,fkaisr,fkaisi,    
     5fkaybg,vbe,wbe,vsoe,wsoe,vise,wise,vbcd(905,2),wbcd(905,2),       
     6vsov(905,2),wsov(905,2),pnloc(2)         
      common/theta/jmax,j1min,j1max,ktheta,mmax,nnnn,
     1thetad(181),sint(181),    
     1cost(181),plmd(181,10),xdwb(181),poldw(181),asymdw(181),
     2tenpol(181,6),plmq(181,10)                     
      common/zf/ktzf(2),kk2,kk3                                         
      common/add/itape                                                  
      common/dbl/dltmas(2),dlpmas(2)                                    
      dimension a(8)          
      character*12 ddsn   
      common/dataset/ddsn                                       
  200 read(35,3050,end=9000) ktout,ins,(numrun(i),i=1,3),ititol,ddsn    
      if(ins.ne.0) numrun(5)=ins                                        
      if(ktout(1).eq.4.or.ktout(1).eq.5.or.ktout(1).eq.8) itape=1       
      if(ktout(1).eq.9) go to 300                                       
 2000 read(35,3000) (a(i),i=1,8)                                        
 3000 format(8f10.5)                                                    
 3050 format(10i1,3i2,i4,12a4,a12)                                      
      go to 4000                                                        
 9000 if(itape.eq.1) close(8)                                         
      if(itape.eq.1) close(8)                                           
      stop
  300 call mix                                                          
      go to 200                                                         
 4000 i=a(1)+0.1                                                        
      pi=i                                                              
      if(i.eq.0) go to 4500                                             
      go to (1,2,3,4,5,6,7,8,9),i                                       
    1 nubchn=a(2)+0.01                                                  
      rmax=a(3)                                                         
      nrmin=a(4)+0.01                                                   
      nrmax=a(5)+0.01                                                   
      elabi=a(6)                                                        
      if(nrmin.eq.0) nrmin=1                                            
      go to 2000                                                        
    2 ich=(a(1)-pi)*10.0+0.01                                           
      trs(1,ich)=a(2)                                                   
      ltr(1,ich)=a(3)+0.01                                              
      trj(1,ich)=a(4)                                                   
      trs(2,ich)=a(5)                                                   
      ltr(2,ich)=a(6)+0.01                                              
      trj(2,ich)=a(7)                                                   
      ttrj=a(8)                                                         
      if(ich.eq.2) ttrj=a(4)                                            
      go  to  2000                                                      
    3 ich=(a(1)-pi)*10.0+0.01                                           
      lmind(ich)=a(2)+0.01                                              
      lmaxd(ich)=a(3)+0.01                                              
      beta(ich)=a(4)                                                    
      pnloc(ich)=a(5)                                                   
      go  to  2000                                                      
    4 ich=(a(1)-pi)*10.0+0.01                                           
      pmas(ich)=a(2)                                                    
      tmas(ich)=a(3)                                                    
      dlpmas(ich)=a(2)                                                  
      dltmas(ich)=a(3)                                                  
      pz(ich)=a(4)                                                      
      tz(ich)=a(5)                                                      
      zz(ich)=pz(ich)*tz(ich)                                           
      pspn(ich)=a(6)                                                    
      tspn(ich)=a(7)                                                    
      qvlue(ich)=a(8)                                                   
      go to 2000                                                        
    5 ich=(a(1)-pi)*10.0+0.01                                           
      vd(ich)=a(2)                                                      
      wd(ich)=a(3)                                                      
      vsod(ich)=a(4)                                                    
      wsod(ich)=a(5)                                                    
      rrd(ich)=a(6)                                                     
      ard(ich)=a(7)                                                     
      rcd(ich)=a(8)                                                     
      go  to  2000                                                      
    6 ich=(a(1)-pi)*10.0+0.01                                           
      rsord(ich)=a(2)                                                   
      asord(ich)=a(3)                                                   
      rsoid(ich)=a(4)                                                   
      asoid(ich)=a(5)                                                   
      go  to  2000                                                      
    7 ich=(a(1)-pi)*10.0+0.01                                           
      csdgd(ich)=a(2)                                                   
      rid(ich)=a(3)                                                     
      aid(ich)=a(4)                                                     
      rgd(ich)=a(5)                                                     
      agd(ich)=a(6)                                                     
      ktop(ich)=3                                                       
      if(abs(a(2)).lt.1.e-5) ktop(ich)=1                                
      if(a(5).gt.0.0.and.a(6).gt.0.0) ktop(ich)=2                       
      go to 2000                                                        
    8 ich=(a(1)-pi)*10.0+0.01                                           
      ktisp(ich)=a(2)+0.01                                              
      visd(ich)=a(3)                                                    
      risrd(ich)=a(4)                                                   
      aisrd(ich)=a(5)                                                   
      wisd(ich)=a(6)                                                    
      risid(ich)=a(7)                                                   
      aisid(ich)=a(8)                                                   
      go to 2000                                                        
    9 jmax=a(2)+0.01                                                    
      dtheta=a(3)                                                       
      thetad(1)=a(4)                                                    
      if(a(3).gt.0.0) go to 100                                         
      read(35,3000) (thetad(jj),jj=1,jmax)                              
      go to 120                                                         
  100 do 110 jj=2,jmax                                                  
  110 thetad(jj)=thetad(jj-1)+dtheta                                    
  120 if(thetad(1)) 130,130,140                                         
  130 j1min=2                                                           
      kt1=1                                                             
      sint(1)=0.0d0                                                     
      cost(1)=1.0d0                                                     
      go to 150                                                         
  140 j1min=1                                                           
      kt1=0                                                             
  150 if(thetad(jmax)-180.0) 170,160,160                                
  160 j1max=jmax-1                                                      
      kt2=3                                                             
      sint(jmax)=0.0d0                                                  
      cost(jmax)=-1.0d0                                                 
      go to 180                                                         
  170 j1max=jmax                                                        
      kt2=1                                                             
  180 ktheta=kt1+kt2                                                    
      do 190 j=j1min,j1max                                              
      th=0.017453293*thetad(j)                                          
      sint(j)=sin(th)                                                   
      cost(j)=cos(th)                                                   
  190 continue                                                          
      go  to  2000                                                      
 4500 do 4550 ich=1,2                                                   
      if(abs(rrd(ich)).lt.1.e-5) rrd(ich)=1.0                           
      if(ard(ich).lt.1.e-5) ard(ich)=1.0                                
      if(abs(rcd(ich)).lt.1.e-5) rcd(ich)=1.0                           
      if(abs(rid(ich)).lt.1.e-5) rid(ich)=1.0                           
      if(aid(ich).lt.1.e-5) aid(ich)=1.0                                
      if(abs(rsord(ich)).lt.1.e-5) rsord(ich)=rrd(ich)                  
      if(asord(ich).lt.1.e-5) asord(ich)=ard(ich)                       
      if(abs(rsoid(ich)).lt.1.e-5) rsoid(ich)=rid(ich)                  
      if(asoid(ich).lt.1.e-5) asoid(ich)=aid(ich)                       
 4550 continue                                                          
      if(tspn(1).lt.0.1) ttrj=tspn(2)                                   
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine inputb                                                 
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/ffc/amp,d02,fts,ta,tb,fnrng,ktff(5),nrange,nafac,nbfac,lsa 
     1,lsb,iread(15),jjj,fread(28),wread(21),ktrl(20),ff(905),ffi(905), 
     2ffone(905),nrangd(2,4),n1xd(2,4),itap(66,10,2),frfac(2,4),frso(2) 
      complex*16 ffone,frfac                                            
      common/zrff/fis,fls,fjs,ffout(10),bb,bsubl,nrn,lss,redwd,zeta,    
     1fma,fmas,fmb,fmbs,sa,fja,sb,fjb,za,zas,zb,zbs,vb(905)             
      dimension  a(8),fmat(15)                                          
 3000 read(35,3050,err=8000) (a(i),i=1,8)                               
 3050 format(8f10.5)                                                    
      i=a(1)                                                            
      if(i) 7000,7000,3100                                              
 3100 pi=i                                                              
      j=(a(1)-pi)*10.0+0.001                                            
      q=j                                                               
      k=((a(1)-pi)*10.0-q)*10.0+0.001                                   
      j1=j+1                                                            
      go to (701,710,720,730,740,750,760,770,780),j1                    
  701 k=k+1                                                             
      go to (702,703),k                                                 
  702 amp=a(2)                                                          
      ktff(4)=a(3)                                                      
      d02=a(4)                                                          
      ktff(5)=a(5)                                                      
      fts=a(6)                                                          
      ta=a(7)                                                           
      tb=a(8)                                                           
      if(abs(d02).lt.0.001) d02=1.018e+04                               
      if(abs(amp).lt.1.e-10) amp=1.0                                    
      go to 3000                                                        
  703 fnrng=a(2)                                                        
      go to 3000                                                        
  710 ktff(1)=1                                                         
      go to (711,712),k                                                 
  711 ktff(2)=a(2)                                                      
      fread(1)=a(3)                                                     
      fread(2)=a(4)                                                     
      fread(3)=a(5)                                                     
      fread(4)=a(6)                                                     
      fread(5)=a(7)                                                     
      go to 3000                                                        
  712 ktff(3)=a(2)                                                      
      fread(6)=a(3)                                                     
      fread(7)=a(4)                                                     
      fread(8)=a(5)                                                     
      fread(9)=a(6)                                                     
      fread(10)=a(7)                                                    
      fread(11)=a(8)                                                    
      if(abs(a(3)).lt.0.001) ktff(3)=0                                  
      go to 3000                                                        
  720 ktff(1)=2                                                         
      go to (721,722),k                                                 
  721 iread(1)=a(2)                                                     
      fread(1)=a(3)                                                     
      fread(2)=a(4)                                                     
      fread(3)=a(5)                                                     
      fread(4)=a(6)                                                     
      go to 3000                                                        
  722 fread(5)=a(2)                                                     
      fread(6)=a(3)                                                     
      fread(7)=a(4)                                                     
      go to 3000                                                        
  730 ktff(1)=3                                                         
      iread(1)=a(2)                                                     
      fread(1)=a(3)                                                     
      go to 3000                                                        
  740 ktff(1)=4                                                         
      go to (741,742,743),k                                             
  741 iread(1)=a(2)                                                     
      iread(2)=a(3)                                                     
      fread(1)=a(4)                                                     
      fread(2)=a(5)                                                     
      fread(3)=a(6)                                                     
      fread(4)=a(7)                                                     
      fread(5)=a(8) 
      frso(1)=-10.d0                                                 
      frso(2)=-10.d0                                                    
      go to 3000                                                        
  742 fread(6)=a(2)                                                     
      fread(7)=a(3)                                                     
      fread(8)=a(4)                                                     
      fread(9)=a(5)                                                     
      fread(23)=a(6)                                                    
      if(fread(7).lt.0.001) fread(7)=1.0                                
      go to 3000               
  743 frso(1)=a(2)                                                     
      frso(2)=a(3)                                                      
      go to 3000                                                        
  750 ktff(1)=5                                                         
      go to(751,752,753,754),k                                          
  751 read(35,3050)(ff(ii),ii=1,nrmax)                                  
      go to 3000                                                        
  752 ktff(3)=20                                                        
      read(35,3050)(ffi(ii),ii=1,nrmax)                                 
      go to 3000                                                        
  753 read(35,3750)(fmat(ii),ii=1,15)                                   
 3750 format(15a4)                                                      
      read  (5,fmat)(ff(ii),ii=1,nrmax)                                 
      go to 3000                                                        
  754 ktff(3)=20                                                        
      read(35,3750)(fmat(ii),ii=1,15)                                   
      read(35,fmat)(ffi(ii),ii=1,nrmax)                                 
      go to 3000                                                        
  760 ktff(1)=6                                                         
      go to (761,762,763,764,765,766),k                                 
  761 iread(3)=a(2)                                                     
      fread(10)=a(3)                                                    
      fread(11)=a(4)                                                    
      go to 3000                                                        
  762 fread(2)=a(2)                                                     
      fread(3)=a(3)                                                     
      go to 3000                                                        
  763 fread(1)=a(2)                                                     
      iread(1)=a(3)                                                     
      iread(4)=a(4)                                                     
      fread(12)=a(5)                                                    
      iread(2)=a(6)                                                     
      fread(4)=a(7)                                                     
      fread(5)=a(8)                                                     
      go to 3000                                                        
  764 fread(6)=a(2)                                                     
      fread(7)=a(3)                                                     
      fread(8)=a(4)                                                     
      fread(9)=a(5)                                                     
      if(fread(7).lt.0.001) fread(7)=1.0                                
      go to 3000                                                        
  765 fread(13)=a(2)                                                    
      iread(5)=a(3)                                                     
      iread(6)=a(4)                                                     
      fread(14)=a(5)                                                    
      iread(10)=a(6)                                                    
      fread(15)=a(7)                                                    
      fread(16)=a(8)                                                    
      go to 3000                                                        
  766 fread(17)=a(2)                                                    
      fread(18)=a(3)                                                    
      fread(19)=a(4)                                                    
      fread(20)=a(5)                                                    
      if(fread(18).lt.0.001) fread(18)=1.0                              
      go to 3000                                                        
  770 ktff(1)=7                                                         
      go to (771,762,763,764,765,766),k                                 
  771 fread(10)=a(2)                                                    
      fread(11)=a(3)                                                    
      go to 3000                                                        
  780 ktff(1)=8                                                         
      go to (781,782),k                                                 
  781 do 7810 nn=1,6                                                    
      nn1=nn+1                                                          
 7810 fread(nn)=a(nn1)                                                  
      go to 3000                                                        
  782 do 7820 nn=7,12                                                   
      nn1=nn-5                                                          
 7820 fread(nn)=a(nn1)                                                  
      go to 3000                                                        
 7000 return                                                            
 8000 write(6,8100)                                                     
 8100 format('0input data error in inputb')                             
      stop
      end                                                               
c     ******************************************************************
      subroutine inte2(ktfg)                                            
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/pot/vd(2),wd(2),vsod(2),wsod(2),rrd(2),ard(2),rcd(2),      
     1csdgd(2),rid(2),aid(2),rsord(2),asord(2),rsoid(2),asoid(2),       
     2ktop(2),rgd(2),agd(2),ktisp(2),risrd(2),aisrd(2),visd(2),risid(2),
     3aisid(2),wisd(2),csd,drhod(12),rhobn,rhoin,rhosor,rhosoi,rhoisr,  
     4rhoisi,rhobc,rhobng,fkayar,fkayai,fkasor,fkasoi,fkaisr,fkaisi,    
     5fkaybg,vbe,wbe,vsoe,wsoe,vise,wise,vbcd(905,2),wbcd(905,2),       
     6vsov(905,2),wsov(905,2),pnloc(2)                                  
      common/coul/z(16),sigld(90,2),fd(90,2),gd(90,2),fpd(90,2),        
     1gpd(90,2)                                                         
      common/damgk2/iromi,irm,liad(90),ljiad(90),lmad(90),ljmad(90),    
     1icd(90),nadli(90,3),cf(905),ci(905,7),xx(905),yy(905),ur,ui,urp,  
     2uip                                                               
      common/zf/ktzf(2),kk2,kk3                                         
c     ***********************************************************       
c     changes for reading wave functions                                
c     ***********************************************************       
      complex*16 wfi(200)
*     complex*16 temp(200)                                     
c     ***********************************************************       
c     changes for spin-dependent wf                                     
c     ***********************************************************       
      complex*16 wfij(200,90,3)                                         
c     ***********************************************************       
      real idata,ifh,ifhp,idatap
*     real*8 tjat                                   
      complex*16 cf,ci                                                  
      dimension   dxx(7),dyy(7)                                         
      dimension idatap(16),ifhp(2),idata(16),ifh(2)                     
c     **********************************************************        
      save indold,wfi,kount,wfij                                        
c     **********************************************************        
      data idata    /4h    ,4h   i,4hncid,4hent ,4h    ,4h    ,4h   e,  
     1 4hxit ,4hinte,4hrmed,4hiate,4h(1) ,4hinte,4hrmed,4hiate,4h(2) /  
      data ifh   /4h  f ,4h  h /                                        
      data tenp35,tenm35/1.0e+35,1.0e-35/                               
c     **********************************************************        
      data indold,kount/77,0/                                           
c     **********************************************************        
      if(kmemo.eq.1972) go to 10                                        
      do 5 k=1,2                                                        
      ifhp(k)=ifh(k)                                                    
    5 idatap(k)=idata(k)                                                
      do 6 k=3,16                                                       
    6 idatap(k)=idata(k)                                                
      kmemo=1972                                                        
   10 flc=lc                                                            
      fic=pspn(ichnl)                                                   
      fl2=flc*(flc+1.0)                                                 
      nv=nint(2.d0*fic+1.d0)
      l=lc+1                                                            
      i=ichnl                                                           
      ich2=i-2                                                          
      fsols=(fjc*(fjc+1.0)-fl2-fic*(fic+1.0))/2.0                       
      drho12=drhod(i+4)                                                 
      drho56=drhod(i+8)                                                 
      rho=drho                                                          
      rho2=rho*rho                                                      
      nrst=1                                                            
      nrnmx=nr3max                                                      
      qr1=0.0                                                           
      qi1=0.0                                                           
      xr1=0.0                                                           
      xi1=0.0                                                           
      if(lc.ne.1) go to 21                                              
      xr1=-rho2/6.0                                                     
      go to 21                                                          
   21 vbs1=vsov(nrst,i)*fsols                                           
      wbs1=wsov(nrst,i)*fsols                                           
      vb1=vbcd(nrst,i)                                                  
      wbc=wbcd(nrst,i)                                                  
      qr2=vb1+vbs1-fl2/rho2                                             
      qi2=wbc+wbs1                                                      
      xr2=rho**l                                                        
      xi2=xr2*rho                                                       
      xx(1)=xr2                                                         
      yy(1)=xi2                                                         
c     **********************************************************        
      if(ktout(i+2).eq.1.or.ktout(i+2).eq.2) go to 788                  
c     **********************************************************        
      do 120 nrn=2,nrnmx                                                
      rho=rho+drho                                                      
      nr=nrn                                                            
      nrm1=nr-1                                                         
      rho2=rho*rho                                                      
      vb1=vbcd(nr,i)+vsov(nr,i)*fsols                                   
      wb1=wbcd(nr,i)+wsov(nr,i)*fsols                                   
      qr3=vb1-fl2/rho2                                                  
      qi3=wb1                                                           
      t1=1.0+drho12*qr3                                                 
      t2=drho12*qi3                                                     
      d1=t2*t2+t1*t1                                                    
      t3=2.0-drho56*qr2                                                 
      t4=-drho56*qi2                                                    
      t5=-(1.0+drho12*qr1)                                              
      t6=-drho12*qi1                                                    
      t8=t3*xr2-t4*xi2+t5*xr1-t6*xi1                                    
      t9=t3*xi2+t4*xr2+t5*xi1+t6*xr1                                    
      xr3=(t8*t1+t9*t2)/d1                                              
      xi3=(-t8*t2+t9*t1)/d1                                             
      if(max(abs(xr3),abs(xi3))-tenp35)  68,68,65                       
   65 continue                                                          
      fac=tenm35                                                        
      xr2=xr2*fac                                                       
      xi2=xi2*fac                                                       
      xr3=xr3*fac                                                       
      xi3=xi3*fac                                                       
      go to 69                                                          
   68 fac=1.0                                                           
   69 xr1=xr2                                                           
      xi1=xi2                                                           
      xr2=xr3                                                           
      xi2=xi3                                                           
      qr1=qr2                                                           
      qi1=qi2                                                           
      qr2=qr3                                                           
      qi2=qi3                                                           
      xx(nr)=xr3                                                        
      yy(nr)=xi3                                                        
      if(fac-0.5)  97,120,120                                           
   97 do 98 nrt=1,nrm1                                                  
      xx(nrt)=xx(nrt)*fac                                               
   98 yy(nrt)=yy(nrt)*fac                                               
  120 continue                                                          
      go to 145                                                         
  145 continue                                                          
      t1=1.0/(60.0*drho)                                                
      t2=9.0*t1                                                         
      t3=5.0*t2                                                         
      do 160 nrt=nr3min,nr3max                                          
      n1=nrt-nr3min+1                                                   
      dxx(n1)=xx(nrt)                                                   
      dyy(n1)=yy(nrt)                                                   
  160 continue                                                          
      ur=dxx(4)                                                         
      ui=dyy(4)                                                         
      urp=t1*(dxx(7)-dxx(1))+t2*(dxx(2)-dxx(6))+t3*(dxx(5)-dxx(3))      
      uip=t1*(dyy(7)-dyy(1))+t2*(dyy(2)-dyy(6))+t3*(dyy(5)-dyy(3))      
      call fnsg2                                                        
  788 if(ktout(i+2).eq.0) go to 7000                                    
      if(ktout(i+2).gt.2) goto 789                                      
c     **********************************************************        
      if(numrun(5).eq.1) then                                           
      indre=16                                                          
      if(indold.ne.indre) kount=0                                       
      indold=16                                                         
      endif                                                             
      if(numrun(5).eq.2) then                                           
      indre=18                                                          
      if(indold.ne.indre) kount=0                                       
      indold=18                                                         
      endif                                                             
      if(numrun(5).eq.3) then                                           
      indre=20                                                          
      if(indold.ne.indre) kount=0                                       
      indold=20                                                         
      endif                                                             
      if(numrun(5).eq.4) then                                           
      indre=22                                                          
      if(indold.ne.indre) kount=0                                       
      indold=22                                                         
      endif                                                             
      if(numrun(5).eq.5) then                                           
      indre=24                                                          
      if(indold.ne.indre) kount=0                                       
      indold=24                                                         
      endif                                                             
      if(numrun(5).eq.6) then                                           
      indre=26                                                          
      if(indold.ne.indre) kount=0                                       
      indold=26                                                         
      endif                                                             
      if(numrun(5).eq.7) then                                           
      indre=28                                                          
      if(indold.ne.indre) kount=0                                       
      indold=28                                                         
      endif                                                             
      if(numrun(5).eq.8) then                                           
      indre=30                                                          
      if(indold.ne.indre) kount=0                                       
      indold=30                                                         
      endif                                                             
      kount=kount+1                                                     
      if(kount.eq.1) then   
      iwff=0                                            
      print *,' reading wavefunctions on channel',indre                 
      rewind indre                                                      
      do 1131 jatl=0,lmaxd(ichnl)                                       
      jatll=jatl+1  
      jlow=1                                                
      if(jatl.eq.0) jlow=nv                                             
      jhigh=nv                                            
      if(jatl.eq.lmax) jhigh=1                                          
      do 1131 jati=jlow,jhigh                                           
      read(indre) wfi    
      iwff=iwff+1                                               
      do 1132 nr=1,nrmax                                                
      wfij(nr,jatll,jati)=wfi(nr+1)                                     
 1132 continue                                                          
 1131 continue 
      print*
      print*,iwff,'  wavefunctions read'                                
      endif  
        if(nv.eq.1) then
        jati=1   
        else                                                      
        jati=int(fjc-flc+2.01) 
        endif                                    
c     cossig=+cos(sigld(lc+1,i))                                        
c     sinsig=+sin(sigld(lc+1,i))                                        
      do 889 nr=1,nrmax                                                 
      xx(nr)=dreal(wfij(nr,lc+1,jati))                                  
      yy(nr)=dimag(wfij(nr,lc+1,jati)) 
c     multiply by coulomb phases if from adia !!      
c     tjat  =xx(nr)*cossig-yy(nr)*sinsig                                  
c     yy(nr)=yy(nr)*cossig+xx(nr)*sinsig                                
c     xx(nr)=tjat                                                       
  889 continue                                                          
  789 if(ktout(i+2).le.1) goto 7000                                     
c     **********************************************************        
      if(iskip1) 510,510,515                                            
  510 iskip1=1                                                          
      go to 520                                                         
  515 if(iskip2-noline+5) 520,520,522                                   
  520 call skip(iskip2,npgs,numrun,ititol)                              
  522 kmin=4*i-3                                                        
      kmax=4*i                                                          
cliu      write(6,525)(idatap(k),k=kmin,kmax),ifhp(ktfg),lc,fjc,nrmax       
  525 format(1h ,2x,4a4,20hchannel distort.wave,a4,2x,                  
     13h(l=,i3,3x,2hj=,f5.1,3x,5hnr=1,,i3,1h))                          
      iskip2=iskip2+2                                                   
cliu      write(6,530)                                                      
  530 format(1h ,2x,9hreal part)                                        
      iskip2=iskip2+2                                                   
      do 550 nr=1,nrmax,5                                               
cliu      write(6,535) xx(nr),xx(nr+1),xx(nr+2),xx(nr+3),xx(nr+4)           
  535 format(1h ,5e21.8)                                                
      iskip2=iskip2+1                                                   
      if(iskip2+1.gt.noline) call skip(iskip2,npgs,numrun,ititol)       
  550 continue                                                          
cliu      write(6,555)                                                      
  555 format(1h ,2x,9himag part)                                        
      iskip2=iskip2+2                                                   
      do 570 nr=1,nrmax,5                                               
cliu      write(6,535) yy(nr),yy(nr+1),yy(nr+2),yy(nr+3),yy(nr+4)           
      iskip2=iskip2+1                                                   
      if(iskip2+1.gt.noline) call skip(iskip2,npgs,numrun,ititol)       
  570 continue                                                          
 7000 return                                                            
      end                                                               
c     ******************************************************************
      subroutine fnsg2                                                  
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/coul/z(16),sigld(90,2),fd(90,2),gd(90,2),fpd(90,2),        
     1gpd(90,2)                                                         
      common/damgk2/iromi,irm,liad(90),ljiad(90),lmad(90),ljmad(90),    
     1icd(90),nadli(90,3),cf(905),ci(905,7),xx(905),yy(905),ur,ui,urp,  
     2uip                                                               
      complex*16 cf,ci                                                  
      dimension fnr(4),fni(4)                                           
      l=lc+1                                                            
      i=ichnl                                                           
      eta=z(i)                                                          
      if(eta-0.0001) 11,12,12                                           
   11 iz1=0                                                             
      go to 13                                                          
   12 iz1=1                                                             
      s6=sigld(l,i)                                                     
      s1=cos(s6)                                                        
      s2=sin(s6)                                                        
   13 f1=fd(l,i)                                                        
      g1=gd(l,i)                                                        
      fp1=fpd(l,i)                                                      
      gp1=gpd(l,i)                                                      
      facu=1.0/urp                                                      
      ur=ur*facu                                                        
      ui=ui*facu                                                        
      uip=uip*facu                                                      
      t1=fp1*ur-f1                                                      
      t2=fp1*ui-f1*uip                                                  
      t3=g1-gp1*ur+t2                                                   
      t4=-g1*uip+gp1*ui+t1                                              
      t5=1.0/(t3*t3+t4*t4)                                              
      t6=(-t2*t4+t1*t3)*t5                                              
      t7=(t1*t4+t2*t3)*t5                                               
      t1=f1+t6*g1-t7*f1                                                 
      t2=t7*g1+t6*f1                                                    
      t3=1.0/(ur*ur+ui*ui)                                              
      t4=(ur*t1+ui*t2)*t3                                               
      t5=(-ui*t1+ur*t2)*t3                                              
      if(iz1-1) 19,18,18                                                
   18 fnr(i)=(t4*s1-t5*s2)*facu                                         
      fni(i)=(t4*s2+t5*s1)*facu                                         
      go to 25                                                          
   19 fnr(i)=t4*facu                                                    
      fni(i)=t5*facu                                                    
   25 do 30 nr=1,nrmax                                                  
      t1=xx(nr)*fnr(i)-yy(nr)*fni(i)                                    
      yy(nr)=yy(nr)*fnr(i)+xx(nr)*fni(i)                                
      xx(nr)=t1                                                         
   30 continue                                                          
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine legp                                                   
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/theta/jmax,j1min,j1max,ktheta,mmax,nnnn,
     1thetad(181),sint(181),    
     1cost(181),plmd(181,10),xdwb(181),poldw(181),asymdw(181),
     2tenpol(181,6),plmq(181,10)
      common/d0out/p(181,10),pp(181,10),t3(181),sq22(181,10),sq33(181)  
      t1=2*lc+1                                                         
      do 150 m1=1,mmax                                                  
      do 150 j=j1min,j1max                                              
      c1=cost(j)                                                        
      if(lc-m1+1) 5000,130,135                                          
  130 if(lc.eq.0) go to 132                                             
      t2=2*lc-1                                                         
      s1=sint(j)                                                        
      t3(j)=t3(j)*t2*s1                                                 
      go to 134                                                         
  132 t3(j)=1.0d0                                                       
  134 t4=t3(j)                                                          
      p(j,m1)=t1*c1*t4                                                  
      go to 140                                                         
  135 t4=p(j,m1)                                                        
      t5=lc-m1+2                                                        
      t6=lc+m1-1                                                        
      p(j,m1)=(t1*c1*t4-t6*pp(j,m1))/t5                                 
  140 pp(j,m1)=t4                                                       
  150 plmd(j,m1)=t4                                                     
      ktrans=ktheta                                                     
      go to (5000,180,190,170),ktrans                                   
  170 plmd(1,1)=1.0                                                     
      plmd(jmax,1)=(-1.0)**lc                                           
      go to 200                                                         
  180 plmd(1,1)=1.0                                                     
      go to 200                                                         
  190 plmd(jmax,1)=(-1.0)**lc                                           
  200 if(mmax-1) 5000,5000,205                                          
  205 go to (5000,220,230,210),ktrans                                   
  210 do 215 m=2,mmax                                                   
      plmd(1,m)=0.0                                                     
      plmd(jmax,m)=0.0                                                  
  215 continue                                                          
      go to 5000                                                        
  220 do 225 m=2,mmax                                                   
      plmd(1,m)=0.0                                                     
  225 continue                                                          
      go to 5000                                                        
  230 do 235 m=2,mmax                                                   
      plmd(jmax,m)=0.0                                                  
  235 continue                                                          
 5000 continue                                                          
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine legq                                                   
      implicit real*8 ( a-h , o-z )                                     
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/theta/jmax,j1min,j1max,ktheta,mmax,nnnn,
     1thetad(181),sint(181),    
     1cost(181),plmd(181,10),xdwb(181),poldw(181),asymdw(181),
     2tenpol(181,6),plmq(181,10)               
      common/d0out/p(181,10),pp(181,10),t3(181),t1(181,10),t2(181)      
      lt=lc+1                                                           
      rl=dfloat(lc)                                                     
      rl1=dfloat(lc-1)                                                  
      rl21=dfloat(2*lc-1)                                               
      rlq=rl+1.0d0                                                      
      if (lt.gt.3) lt=3                                                 
      go to (100,200,300),lt                                            
  100 do 110 m= 1,mmax                                                  
      plmq(1,m)=0.0d0                                                   
      plmq(jmax,m)=0.0d0                                                
  110 continue                                                          
      t2(1)=0.0d0                                                       
      t2(jmax)=0.0d0                                                    
      do 120 k = j1min, j1max                                           
      plmq(k,1)=dlog((1.0d0+cost(k))/(1.0d0-cost(k)))/2.0d0             
      t2(k)=plmq(k,1)                                                   
 120  continue                                                          
      go to 1000                                                        
 200  do 210 m = 1 , mmax                                               
      plmq(1,m)=0.0d0                                                   
      plmq(jmax,m)=0.0d0                                                
      t1(1,m)=0.0d0                                                     
      t1(jmax,m)=0.0d0                                                  
 210  continue                                                          
      do 220 k = j1min , j1max                                          
      plmq(k,1)=cost(k)*t2(k)-1.0d0                                     
      t1(k,1)=plmq(k,1)                                                 
      plmq(k,2)=sint(k)*t2(k)+cost(k)/sint(k)                           
      t1(k,2)=plmq(k,2)                                                 
 220  continue                                                          
      go to 1000                                                        
 300  do 310 m = 1 , mmax                                               
       plmq(1,m)=0.0d0                                                  
       plmq(jmax,m)=0.0d0                                               
 310  continue                                                          
      do 350 k =j1min , j1max                                           
       plmq(k,1)=(rl21/rl)*cost(k)*t1(k,1)-(rl1/rl)*t2(k)               
       do 320 m = 1 , mmax                                              
        if ( lc+1-m .lt. 0) go to 330                                   
        m1=m+1                                                          
        rm=dfloat(m-1)                                                  
        plmq(k,m1)=-((rl-rm)*cost(k)*plmq(k,m)-(rl+rm)*t1(k,m))/sint(k) 
 320   continue                                                         
 330   continue                                                         
       t2(k)=t1(k,1)                                                    
       do 340 im = 1 , mmax                                             
        t1(k,im)=plmq(k,im)                                             
 340   continue                                                         
 350  continue                                                          
 1000 continue                                                          
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine tasktr                                                 
      implicit real*8(a-h,o-z)                                          
      common/bamp/samp(181,10,3,3),samf(181,10,3,3)                     
      complex*16 samp,samf                                              
      do 120 n4=1,3                                                     
      do 120 n3=1,3                                                     
      do 120 n2=1,10                                                    
      do 120 n1=1,181                                                   
      samp(n1,n2,n3,n4)=(0.0,0.0)                                       
  120 samf(n1,n2,n3,n4)=(0.0,0.0)                                       
      call tasizn                                                       
      call xsm                                                          
      call output                                                       
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine null                                                   
      implicit real*8(a-h,o-z)                                          
      common/pot/ra(28),ia(2),rb(4),ib(2),rc(7288)                      
      common/ffc/re(6),ie(25),jjj,drf(28),rf(21),if(20),rg(1810),       
     1 ca(905),ig(1336),cb(8),frso(2)                                   
      common/coul/aseg(916)                                             
      complex*16 ca,cb                                                  
      do 5 i=1,28                                                       
    5 ra(i)=0.0                                                         
      do 10 i=1,2                                                       
      ia(i)=0                                                           
   10 ib(i)=0                                                           
      do 15 i=1,4                                                       
   15 rb(i)=0.0                                                         
      do 20 i=1,4088                                                    
   20 rc(i)=0.0                                                         
      do 25 i=1,6                                                       
   25 re(i)=0.0                                                         
      do 30 i=1,25                                                      
   30 ie(i)=0                                                           
      do 31 i=1,28                                                      
   31 drf(i)=0.0d0                                                      
      do 35 i=1,21                                                      
   35 rf(i)=0.0                                                         
      do 40 i=1,20                                                      
   40 if(i)=0                                                           
      do 45 i=1,916                                                     
   45 aseg(i)=0.0                                                       
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine skip(k,npgs,numrun,ititol)                             
      implicit real*8(a-h,o-z)                                          
      dimension numrun(5),ititol(12)                                    
      npgs=npgs+1                                                       
*     write(6,1510) numrun,ititol,npgs                                  
 1510 format(1h ,/,12h run number=,i2,1h-,i2,1h-,i4,3h  -,i3,3h  -,i5,  
     110x,12a4,16x,5hpage ,i5,//)                                       
      k=0                                                               
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine mix                                                    
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/theta/jmax,j1min,j1max,ktheta,mmax,nnnn,
     1thetad(181),sint(181),    
     1cost(181),plmd(181,10),xdwb(181),poldw(181),asymdw(181),
     2tenpol(181,6),plmq(181,10)  
      common/bamp/samp(181,10,3,3),sampp(181,10,3,3)                    
      complex*16 samp,sampp,faix,cgc,thamp
      complex*16 scamp(181,15,3,15,3),ccgc 
      common/ffd/ititod(12,41),nbdata(41),iiii,xsect(181,2),            
     1faix(41),a(4),xdwbst(181),iname(10)                               
      character inam2*20,infile*20,blanc*1
      character*20 names(10)
      dimension xdwl(15,181)
*     dimension xdwm(15,181)
      dimension inameb(5)                                               
      dimension esp(181),esd(181)                                       
      character*12 ddsn   
      common/clebma/faclog(500)                                         
      common/dataset/ddsn                                       
      common/add/itape                                                  
      data inameb/5*0/                                                  
      blanc=' '
      if(itape.eq.1.or.numrun(5).ne.8) then
      rewind 8                                                          
      jmax=181                                                          
      k=0                                                               
      do 5 j=1,jmax                                                     
    5 xdwbst(j)=0.0                                                     
   10 n=1                                                               
      icount=0                                                          
      do 15 n4=1,3                                                      
      do 15 n3=1,3                                                      
      do 15 n2=1,10                                                     
      do 15 n1=1,jmax                                                   
   15 samp(n1,n2,n3,n4)=(0.0,0.0)                                       
   20 read(35,100,end=150) a,iname                                      
  100 format(4f10.5,10a4)                                               
      nbdata(n)=a(1)+0.01                                               
      ia1=a(1)+0.1                                                      
      if(a(1).le.-0.5) ia1=a(1)-0.1                                     
      faix(n)=cmplx(a(2),a(3))                                          
      if(ia1.eq.-9) go to 80                                            
      if(ia1.le.0) go to 30                                             
      n=n+1                                                             
      go to 20                                                          
   30 nubpro=n-1                                                        
      icoh=ia1                                                          
      nout=a(4)+0.1                                                     
      do 60 n=1,nubpro                                                  
      if(nbdata(n).le.k) rewind 8                                       
      go to 40                                                          
   72 if(icount.eq.1) go to 10                                          
      icount=1                                                          
      rewind8                                                           
   40 read(8,end=72,err=41) k,ititol,thetad,mxsi,mxsf,mxjp,jmax,pspn(1),
     1       pspn(2),fk(1),fk(2),ecmd(1),ecmd(2),isiw,isfw,jpw, 
     2       tspn(1),tspn(2),mshif,ktout(1),trs(1,2),ltr(1,2),
     3       trj(1,2),qvlue(1),qvlue(2),pmas(1),pmas(2),tmas(1),tmas(2) 
   41 read(8,err=42) sampp                                              
   42 continue                                                          
      if(k.ne.nbdata(n)) go to 40                                       
      do 45 i=1,12                                                      
   45 ititod(i,n)=ititol(i)                                             
      do 50 n4=1,3                                                      
      do 50 n3=1,3                                                      
      do 50 n2=1,10                                                     
      do 50 n1=1,jmax                                                   
   50 samp(n1,n2,n3,n4)=samp(n1,n2,n3,n4)+faix(n)*sampp(n1,n2,n3,n4)    
   60 continue                                                          
      if(nout.eq.0) go to 61                                            
      write(9) nout,iname,inameb,thetad,mxsi,mxsf,mxjp,jmax,pspn(1)     
     1        ,pspn(2),fk(1),fk(2),ecmd(1),ecmd(2),isiw,isfw,jpw        
     2        ,tspn(1),tspn(2),mshif,ktout(1),trs(1,2),ltr(1,2)
     3        ,trj(1,2)                                                 
      write(9) samp                                                     
   61 continue                                                          
      call xsm                                                          
      if(icoh.eq.0) go to 93                                            
      do 92 j=1,jmax                                                    
   92 xdwbst(j)=xdwbst(j)+xdwb(j)                                       
   93 continue                                                          
*     write(6,200)                                                      
  200 format(1h ,///,1h ,7x,28hdata no.       mixing factor,/)          
*     write(6,300) (nbdata(n),faix(n),(ititod(i,n),i=1,12),n=1,nubpro)  
  300 format(i12,5x,1h(,f13.5,1h,,f13.5,1h),5x,12a4,/)                  
*     write(6,3000)                                                     
 3000 format(1h ,5x,70h*****( two step cross section, polarization and a
     1nalyzing power )*****,//,2x,5htheta,2x,14hx-sec.(mb/str),2x,10hpol
     2ar(0/0),4x,8ht22(0/0),4x,8ht21(0/0),4x,8ht20(0/0),2x,12h   it11(0/
     30),4x,8ht22(0/0),4x,8ht21(0/0),4x,8ht20(0/0))                     
*     write(6,3100)(thetad(j),xdwb(j),poldw(j),(tenpol(j,i),i=1,3),asymd
*    1w(j),(tenpol(j,i),i=4,6),j=1,jmax)                                
*     write(6,6990)                                                     
      srt2=dsqrt(2.d0)                                                  
      srt3=dsqrt(3.d0)                                                  
      do 6992 jjs=1,jmax                                                
      tenpol(jjs,1)=2.d0/srt3*asymdw(jjs)                               
      tenpol(jjs,2)=tenpol(jjs,4)*srt3-tenpol(jjs,6)/srt2               
      tenpol(jjs,3)=-tenpol(jjs,4)*srt3-tenpol(jjs,6)/srt2              
      tenpol(jjs,4)=srt2*tenpol(jjs,6)                                  
      tenpol(jjs,5)=-srt3*tenpol(jjs,5)                                 
c     ---------------------------------------------------------------   
c     mods for sig+, sig0, sig-                                         
c     ---------------------------------------------------------------   
c     tenpol(jjs,6)=(2.d0*tenpol(jjs,2)+tenpol(jjs,3))/srt3             
      tenpol(jjs,6)=(100.d0+1.5d0*tenpol(jjs,1)+0.5d0*tenpol(jjs,3))    
     +             /300.d0*xdwb(jjs)                                    
c     esp(jjs)=2.d0*(poldw(jjs)-tenpol(jjs,1))                          
      esp(jjs)=(100.d0-tenpol(jjs,3))/300.d0*xdwb(jjs)                  
c     esd(jjs)=3.d0*tenpol(jjs,1)-2.d0*poldw(jjs)                       
      esd(jjs)=(100.d0-1.5d0*tenpol(jjs,1)+0.5d0*tenpol(jjs,3))         
     +             /300.d0*xdwb(jjs)                                    
 6992 continue                                                          
 3100 format(0pf7.2,1pe15.6,0p4f12.4,2x,0p4f12.4)                       
 6990 format(1h ,//,2x,5htheta,3x,10h  ay (0/0),4x,8haxx(0/0),4x,
     18hayy(0/0),
c    14x,8hazz(0/0),4x,8haxz(0/0),4x,8h x2(0/0),4x,8h sp(0/0),4x,
c    18h sd(0/0)) 
     14x,8hazz(0/0),4x,8haxz(0/0),4x,8h sig(+1),4x,8h sig(0) ,4x,
     18h sig(-1)) 
*     write(6,91)(thetad(j),(tenpol(j,i),i=1,6),esp(j),esd(j),j=1,jmax) 
   91 format(0pf7.2,0p5f12.4,2x,3e12.4)                                 
      if(jmax.le.2) go to 103                                           
      sint(1)=sin(0.017453293*thetad(1))                                
      sum=sint(1)*xdwb(1)*(thetad(2)-thetad(1))                         
      jjjj=jmax-1                                                       
      do 101 j=2,jjjj                                                   
      sint(j)=sin(0.017453293*thetad(j))                                
  101 sum=sum+sint(j)*xdwb(j)*(thetad(j+1)-thetad(j-1))*0.5             
      sint(jmax)=sin(0.017453293*thetad(jmax))                          
      sum=sum+sint(jmax)*xdwb(jmax)*(thetad(jmax)-thetad(jmax-1))       
      sum=sum*0.109662                                                  
*     write(6,102) sum                                                  
  102 format(1h ,29h*** integrated cross section=,1pe17.6,8h (mb)***)   
  103 continue                                                          
      do 70 j=1,jmax                                                    
   70 xsect(j,1)=xdwb(j)                                                
      call plm(xsect,1,jmax,thetad)                                     
*     if(nout.ne.0) write(6,52)nout,iname                               
   52 format(1h ,'   *** amp. is written on file 9 *** data number =',  
     1 i3,', ',10a4,' ***')                                             
      if(icoh.ne.-2) go to 10                                           
*     write(6,210)                                                      
  210 format(1h ,'     *****  incoherent sum of cross sections *****')  
*     write(6,3000)                                                     
*     write(6,3200) (thetad(j),xdwbst(j),j=1,jmax)                      
 3200 format(0pf13.2,1pe20.6)                                           
      do 95 j=1,jmax                                                    
   95 xsect(j,1)=xdwbst(j)                                              
      call plm(xsect,1,jmax,thetad)                                     
      do 96 j=1,jmax                                                    
   96 xdwbst(j)=0.0                                                     
      go to 10                                                          
   80 npgs=0                                                            
      return                                                            
  150 stop
      endif
      if(numrun(5).eq.8) then
      print*
      print*,' Reading back stored ampltudes and combining '
      print*
*     -------------------------------------------------------------------
*     Implementation later if needed: i.e. to read Euler angles 
*     -------------------------------------------------------------------   
*     if(mxjp.gt.5) then
*      print*,' Too many m-values for samp definitions -so redimension'
*      stop
*     endif  
*     ------------------------------------------------------------------- 
*     read the Euler angles from front6z generated file
*     open(49,file='rot.'//ddsn,status='unknown') 
*     read(49,*) alpha,betaa,gamaa
*     alpha=3.1415926536d0*alpha
*     betaa=3.1415926536d0*betaa
*     gamaa=3.1415926536d0*gamaa 
*     print*,' Euler angles have been read from file ' 
*     print*,real(alpha),real(betaa),real(gamaa)
*     -------------------------------------------------------------------                                                         
      jmax=181                                                          
      do n5=1,3 
       do n3=1,3                                         
        do n2=1,10                                       
         do n1=1,jmax 
          samp(n1,n2,n3,n5)=(0.0,0.0)
          sampp(n1,n2,n3,n5)=(0.0,0.0)
         enddo
        enddo
       enddo                               
      enddo    
      do n5=1,3 
       do n3=1,3                                         
        do n2=1,15                                       
         do n1=1,jmax 
          do n4=1,15                                         
           scamp(n1,n2,n3,n4,n5)=(0.0,0.0)  
          enddo
         enddo
        enddo
       enddo                               
      enddo                                       
      do j=1,jmax 
       do n1=1,15  
        xdwl(n1,j)=0.0  
*       xdwm(n1,j)=0.0  
       enddo                                                  
      enddo                                                  
      n=1  
*     -------------------------------------------------------------------
*     process and print input summary
*     -------------------------------------------------------------------
  920 read(35,109,end=930) a(1),a(2),inam2                              
  109 format(2f10.5,a20)                                               
      a(3)=0.d0
      ipos=index(inam2,blanc)-1
      names(n)='amp.'//inam2(1:ipos)
      nbdata(n)=a(1)+0.01                                               
      ia1=a(1)+0.1                                                   
      faix(n)=cmplx(a(2),a(3))                                     
      if(ia1.le.0) go to 930                                        
      n=n+1                                                             
      go to 920                                                      
  930 nubpro=n-1                                                       
      write(6,209)                                                      
  209 format(1h ,7x,28hdata no.       mixing factor,/)          
      write(6,309) (nbdata(n),faix(n),names(n),n=1,nubpro)  
  309 format(i12,3x,1h(,f9.5,1h,,f9.5,1h),3x,a20,/)                  
*     -------------------------------------------------------------------
*     read back the amplitudes beta from the earlier calculations
*     -------------------------------------------------------------------
      do 960 n=1,nubpro                                                
      infile=names(n)
      open(8,file=infile,form='unformatted',status='old')
*     -------------------------------------------------------------------
*     read back amplitude information and the amplitude
*     -------------------------------------------------------------------                                   
      read(8,err=941) k,ititol,thetad,mxsi,mxsf,mxjp,jmax,pspn(1)
     1      ,pspn(2),fk(1),fk(2),ecmd(1),ecmd(2),isiw,isfw,jpw  
     2      ,tspn(1),tspn(2),mshif,ktout(1),trs(1,2),ltr(1,2)
     3      ,trj(1,2),qvlue(1),qvlue(2),pmas(1),pmas(2),tmas(1),tmas(2) 
  941 read(8,err=942) sampp 
      print*,' amplitude ',n,' read successfully'
  942 continue   
      close(8) 
*     ------------------------------------------------------------------- 
*     copy the original m>0 matrix elements upward into samp array  
*     anticipating using 2*mxjp elements - the lowest mxjp for m<0
*     ------------------------------------------------------------------- 
      do j=1,jmax
       do msi1=1,mxsi                                      
        do msf1=1,mxsf                                                  
         do mjr1=1,mxjp 
         samp(j,mjr1+mxjp,msf1,msi1)=sampp(j,mjr1,msf1,msi1)       
         enddo
        enddo
       enddo
      enddo    
*     ------------------------------------------------------------------- 
*     retrieve the l and j transfer and spin values 
*     ------------------------------------------------------------------- 
      rjval=trj(1,2) 
      lval=ltr(1,2)
      rsi=pspn(1)
      rsf=pspn(2)
*     ------------------------------------------------------------------- 
*     compute all m<0 elements using the Madison symmetry relation       
*     ------------------------------------------------------------------- 
      do j=1,jmax
       do msi1=1,mxsi   
        rmsi=msi1-rsi-1.d0  
        msip=nint(-rmsi+rsi+1)                                    
        do msf1=1,mxsf   
         rmsf=msf1-rsf-1.d0 
         msfp=nint(-rmsf+rsf+1)                                        
         do mjr1=1,mxjp 
          rmf1=0.5d0*(2*mjr1+mshif)
          mjrp=mxjp+1-mjr1  
          phase=(-1)**nint(lval+rsi-rmsi+rsf-rmsf+rjval+rmf1)
          samp(j,mjrp,msf1,msi1)=phase*sampp(j,mjr1,msfp,msip)       
         enddo
        enddo
       enddo
      enddo
*     ------------------------------------------------------------------- 
*     this completes the calculation of the negative m amplitudes  
*     samp stores the full t-matrix in the Madison coordimates
*     ------------------------------------------------------------------- 
*     cross sections, and also the sum for checking purposes
*     ------------------------------------------------------------------- 
*     fac=0.06332574d0*fk(2)*(2.0*tspn(2)+1.0)/(ecmd(1)*ecmd(2)*fk(1))
*     fac=fac/((2.0*pspn(1)+1.0)*(2.0*tspn(1)+1.0))   
*     ------------------------------------------------------------------- 
*     do mjr1=1,2*mxjp          
*      do j=1,jmax
*       do msi1=1,mxsi                                                 
*        do msf1=1,mxsf                                                 
*         thamp=samp(j,mjr1,msf1,msi1)                                 
*         xdwm(mjr1,j)=xdwm(mjr1,j)+thamp*conjg(thamp)                 
*        enddo  
*       enddo   
*       compute cross section summed over m-substates         
*       xdwm(15,j)=xdwm(15,j)+xdwm(mjr1,j)        
*      enddo 
*     enddo                                                            
*     ------------------------------------------------------------------- 
*     multiply by phase space factor 
*     ------------------------------------------------------------------- 
*     do mjr1=1,15
*      do j=1,jmax                                                    
*       xdwm(mjr1,j)=fac*xdwm(mjr1,j)       
*      enddo
*     enddo 
*     open(30,file='30.'//ddsn,status='unknown') 
*     open(29,file='29.'//ddsn,status='unknown')
*     write(30,3199)(thetad(j),xdwm(15,j),j=1,jmax)    
*     write(29,3198)(thetad(j),(xdwm(k,j),k=1,9),xdwm(15,j),j=1,jmax)   
*     -------------------------------------------------------------------                                                                 
      aa=tspn(1)
      iaa=nint(2*aa)
      bb=tspn(2)
      ibb=nint(2*bb)
      ijval=nint(2*rjval)
*     ------------------------------------------------------------------- 
      print*,' -------------------------------------------------'
      print*,' amplitude ',n,' ',names(n)
      print*,' particle spins',real(rsi),real(rsf)
      print*,' target spins  ',real(aa),real(bb)
      print*,' l, j transfer ',real(ltr(1,2)),real(rjval)
      print*,' -------------------------------------------------'
*     ------------------------------------------------------------------- 
*     reconstruct the scattering amplitude in terms of target spins
*     and sum all amplitudes into scamp(theta,beta,sig2,alpha,sig1)                                 
*     ------------------------------------------------------------------- 
*     also need to take account of statistical factors related to twofnr
*     conventions for the overlaps and the usual factors in fac, below.
*     ------------------------------------------------------------------- 
      statf=sqrt(2*rjval+1.d0)/sqrt(2*bb+1)
*     ------------------------------------------------------------------- 
      do ialph=1,nint(2*aa+1)
      ralph=ialph-aa-1.d0
      iralph=nint(2*ralph)
      do ibet=1,nint(2*bb+1)
      rbet=ibet-bb-1.d0
      irbet=nint(2*rbet)
      do mjr1=1,2*mxjp
      rmf1=mjr1-rjval-1.d0 
      irmf1=nint(2*rmf1)
      cgc=cleb(iaa,iralph,ijval,irmf1,ibb,irbet)
      cgc=cgc*faix(n)*statf
*     ------------------------------------------------------------------- 
*     print*,real(iralph/2.),real(irmf1/2.),real(irbet/2.)
*     print*,cgc
*     print*
*     ------------------------------------------------------------------- 
      do msi1=1,mxsi   
      do msf1=1,mxsf   
      do j=1,jmax
      ccgc=cgc*samp(j,mjr1,msf1,msi1)
      scamp(j,ibet,msf1,ialph,msi1)=ccgc+scamp(j,ibet,msf1,ialph,msi1)
      enddo
      enddo
      enddo
*     ------------------------------------------------------------------- 
      enddo
      enddo
      enddo
  960 continue     
*     ------------------------------------------------------------------- 
*     coherent sum of amplitudes taken - compute the cross sections
*     note the statistical factors for the residue nucleus here. twofnr
*     has a j statistical factor in its definition of the overlap also
*     ------------------------------------------------------------------- 
      fac=0.06332574d0*fk(2)*(2.0*tspn(2)+1.0)/(ecmd(1)*ecmd(2)*fk(1))
      fac=fac/((2.0*pspn(1)+1.0)*(2.0*tspn(1)+1.0))   
*     ------------------------------------------------------------------- 
*     cross sections, and also the sum for checking purposes
      do ibet=1,nint(2*bb+1)
      do j=1,jmax
      do ialph=1,nint(2*aa+1)
      do msi1=1,mxsi
      do msf1=1,mxsf
      thamp=scamp(j,ibet,msf1,ialph,msi1)
      xdwl(ibet,j)=xdwl(ibet,j)+thamp*conjg(thamp)
      enddo
      enddo
      enddo
*     ------------------------------------------------------------------- 
*     compute cross section summed over m-substates         
*     ------------------------------------------------------------------- 
      xdwl(15,j)=xdwl(15,j)+xdwl(ibet,j)
      enddo
      enddo
*     ------------------------------------------------------------------- 
*     multiply by phase space factor 
      do j=1,jmax
      do ibet=1,15
      xdwl(ibet,j)=fac*xdwl(ibet,j)
      enddo
      xdwb(j)=xdwl(15,j)
      enddo
      open(20,file='20.'//ddsn,status='unknown')
      open(19,file='19.'//ddsn,status='unknown')
*     ----------------------------------------------------------------
      write(20,3199)(thetad(j),xdwl(15,j),j=1,jmax)
 3199 format(0pf7.2,1pe15.6)
      write(19,3198)(thetad(j),(xdwl(k,j),k=1,9),xdwl(15,j),j=1,jmax)
 3198 format(f7.2,9e11.4,e12.5)
      print*,' computed substate cross sections in   ','19.'//ddsn
      print*,' computed cm cross section in          ','20.'//ddsn
*--------------------------------------------------------------------- 
*     ddsn is data set name extension from the title line of data set 
*     angles for the kinematics
*---------------------------------------------------------------------  
      do j=1,jmax                                              
      th=0.017453293d0*thetad(j)                                        
      sint(j)=sin(th)                                                   
      cost(j)=cos(th)                                                   
      enddo
      if(thetad(1).lt.1.d-6) then
      sint(1)=0.d0
      cost(1)=1.d0
      endif
      if(abs(thetad(jmax)-180.d0).lt.1.d-6) then
      sint(jmax)=0.d0
      cost(jmax)=-1.d0
      endif
*---------------------------------------------------------------------  
*     output the cross sections, centre of mass and lab, to files
*---------------------------------------------------------------------  
      open(21,file='21.'//ddsn,status='unknown')
      print*,' computed the cm  cross section in     ','21.'//ddsn
      open(22,file='22.'//ddsn,status='unknown')
      print*,' computed the lab cross section in     ','22.'//ddsn
      open(23,file='23.'//ddsn,status='unknown')
      print*,' inverse kinematics (light detected)   ','23.'//ddsn
      open(24,file='24.'//ddsn,status='unknown')
      print*,' inverse kinematics (heavy detected)   ','24.'//ddsn
      print*,'------------------------------------------------------'
*---------------------------------------------------------------------
      write(21,3199)(thetad(j),xdwb(j),j=1,jmax)    
      pie=4.d0*atan(1.d0)
      eden=ecmd(1)/(ecmd(1)+qvlue(2))
*---------------------------------------------------------------------  
*     print*,' pmas(1),pmas(2),ecmd(1),tmas(1),tmas(2),ecm,qvlue(2)'
*     print*,  pmas(1),pmas(2),ecmd(1),tmas(1),tmas(2),ecm,qvlue(2)
*---------------------------------------------------------------------  
*     lab cross section (normal kinematics, light in, light detected)
*     see for example: Schiff p111
      gam1=sqrt(pmas(1)*pmas(2)*eden/(tmas(1)*tmas(2)))
*---------------------------------------------------------------------  
*     lab cross section (inverse kinematics, heavy in, heavy detected)
      gam2=sqrt(tmas(1)*tmas(2)*eden/(pmas(1)*pmas(2)))
*---------------------------------------------------------------------  
*     lab cross section (inverse kinematics, heavy in, light detected)
      gam3=sqrt(tmas(1)*pmas(2)*eden/(pmas(1)*tmas(2)))
*---------------------------------------------------------------------  
*     for case 3 inverse kinematics need to use (180-theta) ...
*     so change sign of cosine in these formula for thet3, fac3 
*     use atan2s to get the right quadrant
*---------------------------------------------------------------------  
      do jjs=1,jmax  
       thet1=atan2(sint(jjs),gam1+cost(jjs))*180.d0/pie  
       fac1=sqrt((1.d0+gam1*gam1+2.d0*gam1*cost(jjs))**3)
       fac1=fac1/abs(1.d0+gam1*cost(jjs))   
*
       thet2=atan2(sint(jjs),gam2+cost(jjs))*180.d0/pie 
       fac2=sqrt((1.d0+gam2*gam2+2.d0*gam2*cost(jjs))**3)
       fac2=fac2/abs(1.d0+gam2*cost(jjs))   
*
       thet3=atan2(sint(jjs),gam3-cost(jjs))*180.d0/pie 
       fac3=sqrt((1.d0+gam3*gam3-2.d0*gam3*cost(jjs))**3)
       fac3=fac3/abs(1.d0-gam3*cost(jjs)) 
*---------------------------------------------------------------------  
       write(22,3189) thet1,fac1*xdwb(jjs)   
       write(23,3189) thet2,fac2*xdwb(jjs)  
       write(24,3189) thet3,fac3*xdwb(jjs)  
*---------------------------------------------------------------------  
      enddo
 3189 format(0pf7.2,1pe15.6)
*---------------------------------------------------------------------  
      sum=0.d0
      if(jmax.le.2) go to 3099                                      
      sum=sint(1)*xdwb(1)*(thetad(2)-thetad(1))                         
      jjjj=jmax-1                                                       
      do 1099 j=2,jjjj                                                 
 1099 sum=sum+sint(j)*xdwb(j)*(thetad(j+1)-thetad(j-1))*0.5             
      sum=sum+sint(jmax)*xdwb(jmax)*(thetad(jmax)-thetad(jmax-1))       
      sum=sum*0.109662                                                  
      write(6,2099) sum                                                 
 2099 format(1h ,29h*** integrated cross section=,1pe17.6,8h (mb)***)   
 3099 continue                                                          
      endif
      stop
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine output                                                 
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/theta/jmax,j1min,j1max,ktheta,mmax,nnnn,
     1thetad(181),sint(181),    
     1cost(181),plmd(181,10),xdwb(181),poldw(181),asymdw(181),
     2tenpol(181,6),plmq(181,10)                                       
      common/bamp/samp(181,10,3,3),samf(181,10,3,3)                     
      common/xsecs/xdwl                                                
      dimension xdwl(15,181)
      dimension esp(181),esd(181)                                       
      complex*16 samp,samf                                              
      character*12 ddsn   
      common/dataset/ddsn                                       
      data ifopen/0/                                                    
      call skip(k,npgs,numrun,ititol)                                   
*     write(6,3000)                                                     
 3000 format(1h ,5x,70h*****( two step cross section, polarization and a
     1nalyzing power )*****,//,2x,5htheta,2x,14hx-sec.(mb/str),2x,10hpol
     2ar(0/0),4x,8ht22(0/0),4x,8ht21(0/0),4x,8ht20(0/0),2x,12h   it11(0/
     30),4x,8ht22(0/0),4x,8ht21(0/0),4x,8ht20(0/0))                     
 6990 format(1h ,//,2x,5htheta,3x,10h  ay (0/0),4x,8haxx(0/0),4x,8hayy(0
     1/0),4x,8hazz(0/0),4x,8haxz(0/0),4x,8h x2(0/0),4x,8h sp(0/0),4x,8h 
     2sd(0/0))                                                          
*     write(6,3100)(thetad(j),xdwb(j),poldw(j),(tenpol(j,i),i=1,3),asymd
*    1w(j),(tenpol(j,i),i=4,6),j=1,jmax)                                
*     write(6,6990)                                                     
      srt2=dsqrt(2.d0)                                                  
      srt3=dsqrt(3.d0)                                                  
      do 6992 jjs=1,jmax                                                
      tenpol(jjs,1)=2.d0/srt3*asymdw(jjs)                               
      tenpol(jjs,2)=tenpol(jjs,4)*srt3-tenpol(jjs,6)/srt2               
      tenpol(jjs,3)=-tenpol(jjs,4)*srt3-tenpol(jjs,6)/srt2              
      tenpol(jjs,4)=srt2*tenpol(jjs,6)                                  
      tenpol(jjs,5)=-srt3*tenpol(jjs,5)                                 
      tenpol(jjs,6)=(2.d0*tenpol(jjs,2)+tenpol(jjs,3))/srt3             
      esp(jjs)=2.d0*(poldw(jjs)-tenpol(jjs,1))                          
 6992 esd(jjs)=3.d0*tenpol(jjs,1)-2.d0*poldw(jjs)                       
 3100 format(0pf7.2,1pe15.6,0p4f12.4,2x,0p4f12.4)                       
*     write(6,91)(thetad(j),(tenpol(j,i),i=1,6),esp(j),esd(j),j=1,jmax) 
   91 format(0pf7.2,0p8f12.4) 
*--------------------------------------------------------------------- 
*     ddsn is data set name extension from the title line of data set 
*---------------------------------------------------------------------  
*     output the cross sections, centre of mass and lab, to files
*---------------------------------------------------------------------  
      open(21,file='21.'//ddsn,status='unknown')
      open(22,file='22.'//ddsn,status='unknown')
      open(23,file='23.'//ddsn,status='unknown')
      open(24,file='24.'//ddsn,status='unknown')
*---------------------------------------------------------------------
      write(21,3199)(thetad(j),xdwb(j),tenpol(j,1),j=1,jmax)    
 3199 format(0pf7.2,1pe15.6,0pf12.4)
      pie=4.d0*atan(1.d0)
      eden=ecmd(1)/(ecmd(1)+qvlue(2))
*---------------------------------------------------------------------  
*     print*,' pmas(1),pmas(2),ecmd(1),tmas(1),tmas(2),ecm,qvlue(2)'
*     print*,  pmas(1),pmas(2),ecmd(1),tmas(1),tmas(2),ecm,qvlue(2)
*---------------------------------------------------------------------  
*     lab cross section (normal kinematics, light in, light detected)
*     see for example: Schiff p111
      gam1=sqrt(pmas(1)*pmas(2)*eden/(tmas(1)*tmas(2)))
*---------------------------------------------------------------------  
*     lab cross section (inverse kinematics, heavy in, heavy detected)
      gam2=sqrt(tmas(1)*tmas(2)*eden/(pmas(1)*pmas(2)))
*---------------------------------------------------------------------  
*     lab cross section (inverse kinematics, heavy in, light detected)
      gam3=sqrt(tmas(1)*pmas(2)*eden/(pmas(1)*tmas(2)))
*---------------------------------------------------------------------  
*     for case 3 inverse kinematics need to use (180-theta) ...
*     so change sign of cosine in these formula for thet3, fac3 
*     use atan2s to get the right quadrant
*---------------------------------------------------------------------  
      do jjs=1,jmax  
       thet1=atan2(sint(jjs),gam1+cost(jjs))*180.d0/pie  
       fac1=sqrt((1.d0+gam1*gam1+2.d0*gam1*cost(jjs))**3)
       fac1=fac1/abs(1.d0+gam1*cost(jjs))   
*
       thet2=atan2(sint(jjs),gam2+cost(jjs))*180.d0/pie 
       fac2=sqrt((1.d0+gam2*gam2+2.d0*gam2*cost(jjs))**3)
       fac2=fac2/abs(1.d0+gam2*cost(jjs))   
*
       thet3=atan2(sint(jjs),gam3-cost(jjs))*180.d0/pie 
       fac3=sqrt((1.d0+gam3*gam3-2.d0*gam3*cost(jjs))**3)
       fac3=fac3/abs(1.d0-gam3*cost(jjs)) 
*---------------------------------------------------------------------  
*      original lines (jat)
*---------------------------------------------------------------------  
*      thet1=atan(sint(jjs)/(gam1+cost(jjs)))*180.d0/pie    
*      thet2=atan(sint(jjs)/(gam2+cost(jjs)))*180.d0/pie 
*      thet3=atan(sint(jjs)/(gam3+cost(jjs)))*180.d0/pie 
*      fac1=sqrt((1.d0+gam1*gam1+2.d0*gam1*cost(jjs))**3)
*      fac1=fac1/abs(1.d0+gam1*cost(jjs))   
*      fac2=sqrt((1.d0+gam2*gam2+2.d0*gam2*cost(jjs))**3)
*      fac2=fac2/abs(1.d0+gam2*cost(jjs))   
*      fac3=sqrt((1.d0+gam3*gam3+2.d0*gam3*cost(jjs))**3)
*      fac3=fac3/abs(1.d0+gam3*cost(jjs))   
*---------------------------------------------------------------------  
       write(22,3189) thet1,fac1*xdwb(jjs)   
       write(23,3189) thet2,fac2*xdwb(jjs)  
       write(24,3189) thet3,fac3*xdwb(jjs)  
*---------------------------------------------------------------------  
      enddo
 3189 format(0pf7.2,1pe15.6)
*---------------------------------------------------------------------  
      if(jmax.le.2) go to 30                                            
      sum=sint(1)*xdwb(1)*(thetad(2)-thetad(1))                         
      jjjj=jmax-1                                                       
      do 10 j=2,jjjj                                                    
   10 sum=sum+sint(j)*xdwb(j)*(thetad(j+1)-thetad(j-1))*0.5             
      sum=sum+sint(jmax)*xdwb(jmax)*(thetad(jmax)-thetad(jmax-1))       
      sum=sum*0.109662                                                  
      write(6,20) sum                                                   
   20 format(1h ,29h*** integrated cross section=,1pe17.6,8h (mb)***)   
   30 continue                                                          
      if(ktout(7).ne.0) go to 80                                        
      call plm(xdwb,1,jmax,thetad)                                      
*---------------------------------------------------------------------  
*  80 if(ktout(1).eq.0.or.ktout(1).eq.8) go to 7000                                      
*---------------------------------------------------------------------  
   80 if(ktout(1).eq.0) go to 7000                                      
      jjmx=(jmax+1)/2                                                   
      ko4=0                                                             
      ko9=0                                                             
*---------------------------------------------------------------------  
      if(ktout(1).eq.5.or.ktout(1).eq.8) go to 160        
*---------------------------------------------------------------------  
      if(ktout(1).eq.1.or.ktout(1).eq.3.or.ktout(1).eq.4) ko4=1         
      if(ktout(1).eq.2.or.ktout(1).eq.3) ko9=1                          
      if(ko4.eq.0) go to 90                                             
      call skip(k,npgs,numrun,ititol)                                   
   90 if(ko4.eq.1) continue
*     write(6,3310)                                        
 3310 format(1h ,29hbeta(theta,mia,mib,mj) values)                      
      k=k+2                                                             
      pma1=mxsi                                                         
      pma1=(1.0-pma1)/2.0                                               
      pmb1=mxsf                                                         
      pmb1=(1.0-pmb1)/2.0                                               
      do 155 mia1=1,mxsi                                                
      mia=2*mia1-isiw-2                                                 
      fmia1=mia1-1                                                      
      pma2=pma1+fmia1                                                   
      do 155 mib1=1,mxsf                                                
      mib=2*mib1-isfw-2                                                 
      fmib1=mib1-1                                                      
      pmb2=pmb1+fmib1                                                   
      do 155 mjr1=1,mxjp                                                
      mj=2*mjr1+mshif                                                   
      if(mj) 155,95,95                                                  
   95 pmj=float(mj)*0.5                                                 
      if(ko4.eq.1) continue
*     write(6,3320)pma2,pmb2,pmj                           
 3320 format(1h ,5x,4hmia=,f4.1,3x,4hmib=,f4.1,3x,3hmj=,f4.1/10x,
     + 5htheta,13x
     1,7hre.samp,13x,7him.samp,15x,5htheta,13x,7hre.samp,13x,7him.samp) 
      if(ko9.eq.1) write(7,3325) pma2,pmb2,pmj                          
 3325 format(f6.1,f10.1,f9.1)                                           
      k=k+3                                                             
      if(k-noline) 110,105,105                                          
  105 call skip(k,npgs,numrun,ititol)                                   
  110 do 150 j=1,jjmx                                                   
      jj=j+jjmx                                                         
  130 k=k+1                                                             
      if(k-noline) 140,135,135                                          
  135 call skip (k,npgs,numrun,ititol)                                  
  140 if(ko4.eq.1) continue
*     write(6,3330) thetad(j),samp(j,mjr1,mib1,mia1),      
*    1                           thetad(jj),samp(jj,mjr1,mib1,mia1)     
  145 if(ko9.eq.1) write(7,3335) thetad(j),samp(j,mjr1,mib1,mia1),      
     1                           thetad(jj),samp(jj,mjr1,mib1,mia1)     
 3335 format(2(f8.2,2x,2e15.7))                                         
 3330 format(2(f15.2,5x,2e20.8))                                        
  150 continue                                                          
  155 continue                                                          
  160 continue                                                          
      if(ktout(1).ne.4.and.ktout(1).ne.5.and.ktout(1).ne.8) go to 7000
*---------------------------------------------------------------------
c     if(ifopen.eq.0) open(8,status='scratch',form='unformatted',       
c    +   access='sequential')                                           
c     ifopen=5                                                          
*---------------------------------------------------------------------
*     write the amplitude to a file
*---------------------------------------------------------------------
      open(8,file='amp.'//ddsn,form='unformatted',status='unknown')
      write(8) numrun(5),ititol,thetad,mxsi,mxsf,mxjp,jmax,pspn(1),     
     1       pspn(2),fk(1),fk(2),ecmd(1),ecmd(2),isiw,isfw,jpw,       
     2       tspn(1),tspn(2),mshif,ktout(1),trs(1,2),ltr(1,2),
     3       trj(1,2),qvlue(1),qvlue(2),pmas(1),pmas(2),tmas(1),tmas(2) 
      write(8) samp                                                     
      write(6,*)                                                     
      write(6,3345)                                                     
 3345 format(1h ,'*** samp has been written on disk ***')            
*---------------------------------------------------------------------
*     numfar=numrun(5)+10                                               
*     write(8) numfar,ititol,thetad,mxsi,mxsf,mxjp,jmax,pspn(1),        
*    1         pspn(2),fk(1),fk(2),ecmd(1),ecmd(2),isiw,isfw,jpw,       
*    2         tspn(1),tspn(2),mshif                                    
*     write(8) samf                                                     
*     write(6,3345)                                                     
*3345 format(1h ,5x,'*** samn and samf have been written on tape ***')  
*---------------------------------------------------------------------
 7000 return                                                            
      end                                                               
c     ******************************************************************
      subroutine ovlapa(ovlpi,mjtlf,mjtljf,mjtli,mjtlji,iic,ktzf3)      
c     subroutine ovlapa(ovlpi,lf,ljf,li,lji,ic,ktzf3)                   
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/pot/vd(2),wd(2),vsod(2),wsod(2),rrd(2),ard(2),rcd(2),      
     1csdgd(2),rid(2),aid(2),rsord(2),asord(2),rsoid(2),asoid(2),       
     2ktop(2),rgd(2),agd(2),ktisp(2),risrd(2),aisrd(2),visd(2),risid(2),
     3aisid(2),wisd(2),csd,drhod(12),rhobn,rhoin,rhosor,rhosoi,rhoisr,  
     4rhoisi,rhobc,rhobng,fkayar,fkayai,fkasor,fkasoi,fkaisr,fkaisi,    
     5fkaybg,vbe,wbe,vsoe,wsoe,vise,wise,vbcd(905,2),wbcd(905,2),       
     6vsov(905,2),wsov(905,2),pnloc(2)                                  
      common/coul/z(16),sigld(90,2),fd(90,2),gd(90,2),fpd(90,2),        
     1gpd(90,2)                                                         
      common/damgk2/iromi,irm,liad(90),ljiad(90),lmad(90),ljmad(90),    
     1icd(90),nadli(90,3),cf(905),ci(905,7),xx(905),yy(905),ur,ui,urp,  
     2uip                                                               
      complex*16 cf,ci,ffone,frfac,zcf,ovlpi                            
      common/ffc/amp,d02,fts,ta,tb,fnrng,ktff(5),nrange,nafac,nbfac,lsa 
     1,lsb,iread(15),jjj,fread(28),wread(21),ktrl(20),ff(905),ffi(905), 
     2ffone(905),nrangd(2,4),n1xd(2,4),itap(66,10,2),frfac(2,4),frso(2) 
      common/mjt/ zcf(905)                                              
      data lfmem,ljfmem,limem/3*0/                                      
      lf=mjtlf                                                          
      ljf=mjtljf                                                        
      li=mjtli                                                          
      lji=mjtlji                                                        
      ic=iic                                                            
      ljfm=ljf                                                          
      ljim=lji                                                          
      if(abs(vsod(1)).lt.0.0001)  lji=1                                 
      if(abs(vsod(2)).lt.0.0001)  ljf=1                                 
      ljfchk=ljf+100                                                    
      lfchk=lf+100                                                      
      iflagf=0                                                          
      if((ljfchk.eq.ljfmem).and.(lfchk.eq.lfmem))   go to 70            
      ljfmem=ljfchk                                                     
      lfmem=lfchk                                                       
      lc=lf-1                                                           
      ichnl=2                                                           
      drho=drhod(2)                                                     
      fkay=fk(2)                                                        
      fjl=ljf+lf-2                                                      
      fjc=fjl-pspn(2)                                                   
      call  inte2(1)                                                    
      iflagf=1                                                          
      do  60  nr=    1,nrmax                                            
   60 cf(nr)=cmplx(xx(nr),yy(nr))                                       
   70 continue                                                          
      ni=nadli(li,lji)                                                  
      if(ni)   100,100,200                                              
  100 if(iromi.gt.7) iromi=1                                            
      ichnl=1                                                           
      lc=li-1                                                           
      drho=drhod(1)                                                     
      fkay=fk(1)                                                        
      fjl=lji+li-2                                                      
      fjc=fjl-pspn(1)                                                   
      call  inte2(1)                                                    
      do  110  nr=    1,nrmax                                           
  110 ci(nr,iromi)=cmplx(xx(nr),yy(nr))                                 
      lili=liad(iromi)                                                  
      if(lili.eq.0)  go to 112                                          
      ljilji=ljiad(iromi)                                               
      nadli(lili,ljilji)=0                                              
  112 nadli(li,lji)=iromi                                               
      liad(iromi)=li                                                    
      ljiad(iromi)=lji                                                  
      ni=iromi                                                          
      iromi=iromi+1                                                     
  200 continue                                                          
      lichk=li+100                                                      
      if(iflagf.eq.0.and.lichk.eq.limem) go to 555                      
      if(iflagf.eq.0.and.ktzf3.eq.0) go to 555                          
      limem=lichk                                                       
      tra=drd(1)*float(nrmin)                                           
      do 552 nra=nrmin,nrmax                                            
      ntra=1-nrmin+nra                                                  
      zcf(ntra)=cf(nra)*ffone(nra)/tra                                  
  552 tra=tra+drd(1)                                                    
  555 continue                                                          
      ntra=1                                                            
      tra=drd(1)*float(nrmin)                                           
      ovlpi=(0.,0.)                                                     
      do 560 nra=nrmin,nrmax                                            
      k=nra                                                             
      if(nrmin.ne.1) k=nra-nrmin                                        
      coef=0.6666667*float(1+k-2*(k/2))                                 
      if(k.eq.0.or.nra.eq.nrmax)   coef=0.3333333                       
      ovlpi=ovlpi+zcf(ntra)*ci(nra,ni)*tra*coef                         
c     write(1,7777)nra,zcf(ntra)*tra,ci(nra,ni)                         
c7777 format(1h ,i4,3x,2d12.4,3x,2d12.4)                                
      ntra=ntra+1                                                       
  560 tra=tra+drd(1)                                                    
      ovlpi=ovlpi*drd(1)                                                
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine   pgen2                                                
      implicit real*8(a-h,o-z)                                          
c     modified pgen2 to read channel potentials from cards, as well     
c     as using the conventional functional forms                        
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/pot/vd(2),wd(2),vsod(2),wsod(2),rrd(2),ard(2),rcd(2),      
     1csdgd(2),rid(2),aid(2),rsord(2),asord(2),rsoid(2),asoid(2),       
     2ktop(2),rgd(2),agd(2),ktisp(2),risrd(2),aisrd(2),visd(2),risid(2),
     3aisid(2),wisd(2),csd,drhod(12),rhobn,rhoin,rhosor,rhosoi,rhoisr,  
     4rhoisi,rhobc,rhobng,fkayar,fkayai,fkasor,fkasoi,fkaisr,fkaisi,    
     5fkaybg,vbe,wbe,vsoe,wsoe,vise,wise,vbcd(905,2),wbcd(905,2),       
     6vsov(905,2),wsov(905,2),pnloc(2)                                  
      dimension fmat(15)                                                
      data expmax/174.0/                                                
 6000 format (15a4)                                                     
 8000 format('     read in ',i3,' points for chann',i2,' in steps of',  
     1 f7.4,' with format ',2a4,' for real')                            
 8100 format('     read in ',i3,' points for chann',i2,' in steps of',  
     1 f7.4,' with format ',2a4,' for imag')                            
 8200 format('     read in ',i3,' points for chann',i2,' in steps of',  
     1 f7.4,' with format ',2a4,' for reso')                            
 8300 format('     read in ',i3,' points for chann',i2,' in steps of',  
     1 f7.4,' with format ',2a4,' for imso')                            
      grho=drho/fkay                                                    
      ig1=idint(rhobn/50.d0)                                            
      ig2=idint(rhoin/50.d0)                                            
      ig3=idint(rhosor/50.d0)                                           
      ig4=idint(rhosoi/50.d0)                                           
      if(ig1.eq.0) go to 1032 
*     write(6,*)                                          
      write(6,*)                                          
      read(35,6000)(fmat(ii),ii=1,15)                                   
      write(6,8000)nr3max,ichnl,grho,(fmat(ii),ii=1,2 )                 
      read(35,fmat)(vbcd(nr,ichnl),nr=1,nr3max)                         
 1032 if(ig2.eq.0) go to 1033                                           
      read(35,6000)(fmat(ii),ii=1,15)                                   
      read(35,fmat)(wbcd(nr,ichnl),nr=1,nr3max)                         
      write(6,8100) nr3max,ichnl,grho,(fmat(ii),ii=1,2)                 
 1033 if(ig3.eq.0) go to 1034                                           
      read(35,6000)(fmat(ii),ii=1,15)                                   
      read(35,fmat)(vsov(nr,ichnl),nr=1,nr3max)                         
      write(6,8200)nr3max,ichnl,grho,(fmat(ii),ii=1,2)                  
      write(6,*)                                          
 1034 if(ig4.eq.0) go to 1035                                           
      read(35,6000)(fmat(ii),ii=1,15)                                   
      read(35,fmat)(wsov(nr,ichnl),nr=1,nr3max)                         
      write(6,8300)nr3max,ichnl,grho,(fmat(ii),ii=1,2)                  
      write(6,*)                                          
 1035 if(ig3.ne.0) go to 610                                            
      do 600 nr=1,nr3max                                                
      vsov(nr,ichnl)=0.0                                                
  600 continue                                                          
      go to 630                                                         
  610 do 620 nr=1,nr3max                                                
      vsov(nr,ichnl)=-vsov(nr,ichnl)*vsoe                               
  620 continue                                                          
  630 if(ig4.ne.0) go to 650                                            
      do 640 nr=1,nr3max                                                
      wsov(nr,ichnl)=0.0                                                
  640 continue                                                          
      go to 670                                                         
  650 do 660 nr=1,nr3max                                                
      wsov(nr,ichnl)=-wsov(nr,ichnl)*wsoe                               
  660 continue                                                          
  670 continue                                                          
      rhobc2=rhobc*rhobc                                                
      etabc=eta/rhobc                                                   
      fkay4=fkay*fkay*4.0                                               
      rho=drho                                                          
      rho2=rho*rho                                                      
      ffisr=0.0                                                         
      ffisi=0.0                                                         
      dffisr=0.0                                                        
      dffisi=0.0                                                        
      ktisr=ktisp(ichnl)/10                                             
      ktisi=ktisp(ichnl)-ktisr*10                                       
      do  300  nr=1,nr3max                                              
      if(ig1.ne.0) go to 680                                            
      xr=(rho-rhobn)/fkayar                                             
      if(xr.gt.expmax) xr=expmax                                        
      exr=exp(xr)                                                       
      sr1=1.0/(1.0+exr)                                                 
      ffcr=vbe*sr1                                                      
      go to 690                                                         
  680 ffcr=-vbcd(nr,ichnl)*vbe                                          
  690 if(ig2.ne.0) go to 700                                            
      xi=(rho-rhoin)/fkayai                                             
      if(xi.gt.expmax) xi=expmax                                        
      exi=exp(xi)                                                       
      si1=1.0/(1.0+exi)                                                 
      if(ktop(ichnl)-2)  25,30,40                                       
   25 ffci=wbe*si1                                                      
      go  to  50                                                        
   30 sig=exp(-((rho-rhobng)/fkaybg)**2)                                
      ffci=wbe*((1.0-csd)*si1+csd*sig)                                  
      go  to  50                                                        
   40 si2=si1*si1*exi*4.0                                               
      ffci=wbe*((1.0-csd)*si1+csd*si2)                                  
      go to 50                                                          
  700 ffci=-wbcd(nr,ichnl)*wbe                                          
   50 zzz=zz(ichnl)                                                     
      if(zzz-1.0)  60,55,55                                             
   55 if(rho-rhobc)  57,57,59                                           
   57 ffcl=etabc*(3.0-rho2/rhobc2)                                      
      go  to  62                                                        
   59 ffcl=2.0*eta/rho                                                  
      go  to  62                                                        
   60 ffcl=0.0                                                          
   62 if(ig3.ne.0) go to 710                                            
      if(abs(vsoe).lt.0.01)  go  to  65                                 
      xsor=(rho-rhosor)/fkasor                                          
      if(xsor.gt.expmax) xsor=expmax                                    
      esor=exp(xsor)                                                    
      ssor1=1.0/(1.0+esor)                                              
      ssor2=ssor1*ssor1*esor                                            
      vsov(nr,ichnl)=ssor2*fkay4/(fkasor*rho)*vsoe                      
  710 if(ig4.ne.0) go to 720                                            
   65 if(abs(wsoe).lt.0.01)   go  to  41                                
      xsoi=(rho-rhosoi)/fkasoi                                          
      if(xsoi.gt.expmax) xsoi=expmax                                    
      esoi=exp(xsoi)                                                    
      ssoi1=1.0/(1.0+esoi)                                              
      ssoi2=ssoi1*ssoi1*esoi                                            
      wsov(nr,ichnl)=ssoi2*fkay4/(fkasoi*rho)*wsoe                      
  720 continue                                                          
   41 if(ktisr.eq.0)  go to 45                                          
      xisr=(rho-rhoisr)/fkaisr                                          
      if(xisr.gt.expmax) xisr=expmax                                    
      eisr=exp(xisr)                                                    
      sisr1=1.0/(1.0+eisr)                                              
      go to (42,43),ktisr                                               
   42 ffisr=vise*sisr1                                                  
      go to 45                                                          
   43 sisr2=sisr1*sisr1*eisr                                            
      ffisr=vise*sisr2*rhoisr/fkaisr                                    
   45 if(ktisi.eq.0)  go to 80                                          
      xisi=(rho-rhoisi)/fkaisi                                          
      if(xisi.gt.expmax) xisi=expmax                                    
      eisi=exp(xisi)                                                    
      sisi1=1.0/(1.0+eisi)                                              
      go to (46,47),ktisi                                               
   46 ffisi=wise*sisi1                                                  
      go to 80                                                          
   47 sisi2=sisi1*sisi1*eisi                                            
      ffisi=wise*sisi2*rhoisi/fkaisi                                    
   80 vbcd(nr,ichnl)=1.0-ffcl+ffcr+ffisr*tautau(ichnl)                  
      wbcd(nr,ichnl)=ffci+ffisi*tautau(ichnl)                           
      rho=rho+drho                                                      
      rho2=rho*rho                                                      
  300 continue                                                          
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine sigma(lmax,eta,sigma0,sigmad)                          
      implicit real*8(a-h,o-z)                                          
      dimension sigmad(90)                                              
      if(eta.ge.10.) go to 20                                           
      eta2=eta*eta                                                      
      t1=eta2+16.0                                                      
      t2=t1*t1                                                          
      t4=t2*t2                                                          
      t6=eta/6.0                                                        
      sigma0=-(eta/(12.0*t1))*(1.0+(eta2-48.0)/(30.0*t2)+(eta2**2-160.0*
     1eta2+1280.0)/(t4*105.0))-eta+3.0*t6*log(t1)+3.5*atan (1.5*t6)-at  
     2an (eta)-atan (3.0*t6)-atan (2.0*t6)                              
      go to 25                                                          
   20 t1=1.0/eta                                                        
      t2=t1*t1                                                          
      t3=t1*t2                                                          
      t5=t3*t2                                                          
      t7=t5*t2                                                          
      t9=t7*t2                                                          
      sigma0=0.7853981634d0+eta*log(eta)-eta                            
     1      -(0.08333333333d0*t1+0.00277777777d0*t3                     
     2       +0.00079365079d0*t5+0.00059523810d0*t7                     
     3       +0.00084175084d0*t9)                                       
   25 mn=sigma0/6.2831853072d0                                          
      sigma0=sigma0-6.2831853072d0*mn                                   
      sigmad(1)=sigma0                                                  
      do 30 l=2,lmax                                                    
      fl1=l-1                                                           
      sigmad(l)=sigmad(l-1)+atan (eta/fl1)                              
   30 continue                                                          
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine   tasizn                                               
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/pot/vd(2),wd(2),vsod(2),wsod(2),rrd(2),ard(2),rcd(2),      
     1csdgd(2),rid(2),aid(2),rsord(2),asord(2),rsoid(2),asoid(2),       
     2ktop(2),rgd(2),agd(2),ktisp(2),risrd(2),aisrd(2),visd(2),risid(2),
     3aisid(2),wisd(2),csd,drhod(12),rhobn,rhoin,rhosor,rhosoi,rhoisr,  
     4rhoisi,rhobc,rhobng,fkayar,fkayai,fkasor,fkasoi,fkaisr,fkaisi,    
     5fkaybg,vbe,wbe,vsoe,wsoe,vise,wise,vbcd(905,2),wbcd(905,2),       
     6vsov(905,2),wsov(905,2),pnloc(2)                                  
      common/theta/jmax,j1min,j1max,ktheta,mmax,nnnn,
     1thetad(181),sint(181),    
     1cost(181),plmd(181,10),xdwb(181),poldw(181),asymdw(181),
     2tenpol(181,6),plmq(181,10)                     
      common/bamp/samp(181,10,3,3),samf(181,10,3,3)                     
      complex*16 samp,cf,ci,ffone,frfac,samf                            
      common/damgk2/iromi,irm,liad(90),ljiad(90),lmad(90),ljmad(90),    
     1icd(90),nadli(90,3),cf(905),ci(905,7),xx(905),yy(905),ur,ui,urp,  
     2uip                                                               
      common/clebma/faclog(500)                                         
      common/ffc/amp,d02,fts,ta,tb,fnrng,ktff(5),nrange,nafac,nbfac,lsa 
     1,lsb,iread(15),jjj,fread(28),wread(21),ktrl(20),ff(905),ffi(905), 
     2ffone(905),nrangd(2,4),n1xd(2,4),itap(66,10,2),frfac(2,4),frso(2) 
      common/zf/ktzf(2),kk2,kk3                                         
      complex*16 fac7,fac8,fac9,betn,betf,ovlpi,fac3,phase              
      iromi=1                                                           
      pie=4.d0*datan(1.d0)                                              
      irm=1                                                             
      do 10 n1=1,90                                                     
      liad(n1)=0                                                        
      ljiad(n1)=0                                                       
      lmad(n1)=0                                                        
      ljmad(n1)=0                                                       
      icd(n1)=0                                                         
      do 10 n2=1,3                                                      
   10 nadli(n1,n2)=0                                                    
      ktzf3=ktzf(2)                                                     
      if(nubchn.eq.0) ktzf3=ktzf(1)                                     
      ktzf2=ktzf(1)                                                     
      ichnl=1                                                           
      isiw=pspn(1)*2.0+0.01                                             
      isfw=pspn(2)*2.0+0.01                                             
      jtw=tspn(1)*2.0+0.01                                              
      jrw=tspn(2)*2.0+0.01                                              
      jpw=ttrj*2.0+0.01                                                 
      mxsi=isiw+1                                                       
      mxsf=isfw+1                                                       
      mxjp=jpw/2+1                                                      
      mshif=jpw-2*(jpw/2)-2                                             
      limin=lmind(1)+1                                                  
      lfmin=lmind(2)+1                                                  
      limax=lmaxd(1)+1                                                  
      lfmax=lmaxd(2)+1                                                  
      if(lfmin.eq.1)  go to 30                                          
      lfcrmx=lfmin-1                                                    
      do 25 lfcr=1,lfcrmx                                               
      lc=lfcr-1                                                         
      mmax=min0((jpw+isiw+isfw)/2+1,lc+1)                               
      call legq                                                         
   25 call legp                                                         
   30 do 820 lf=lfmin,lfmax                                             
      lfw=(lf-1)*2                                                      
      flfw1=lfw+1                                                       
      fac0=sqrt(flfw1)*12.56637062                                      
      lc=lf-1                                                           
      mmax=min0((jpw+isiw+isfw)/2+1,lc+1)                               
      call legq                                                         
      call legp                                                         
      jfwmin=iabs(lfw-isfw)                                             
      do 810 jf=1,mxsf                                                  
      jfw=(jf+lf-2)*2-isfw                                              
      if(jfw.lt.jfwmin)  go to 810                                      
      fac1=fac0*sqrt(float(jfw+1))                                      
      nubch2=nubchn+2                                                   
      do 800 nch=2,nubch2                                               
      if(abs(beta(nch)).lt.1.e-10) go to 800                            
      ic=nch-2                                                          
      is1w=trs(1,2)*2.0+0.01                                            
      ltr1=ltr(1,2)                                                     
      ltr1w=ltr1*2                                                      
      jtr1w=trj(1,2)*2.0+0.01                                           
      if(ktzf3.eq.0) n1xd(1,2)=0                                        
      limn=iabs(lf-1-ltr1)+1+n1xd(1,2)                                  
      limax1=lf+ltr1                                                    
      limx=min0(limax1,limax)                                           
      if(limx.lt.limn) go to 800                                        
      t1=(is1w+1)*(ltr1w+1)                                             
      fac2=fac1*sqrt(t1)                                                
      do  710  li=limn,limx,2                                           
      if(li.lt.limin) go to 710                                         
      liw=(li-1)*2                                                      
      t1=liw+1                                                          
      fac3=fac2*sqrt(flfw1*t1)*cleb(lfw,0,ltr1w,0,liw,0)*0.282094792    
      phase=(0.d0,1.d0)**((liw-lfw-2*ltr(1,nch))/2)                     
      fac3=fac3*phase                                                   
      jimin=max0(iabs(jfw-jtr1w),iabs(liw-isiw))                        
      jimax=jfw+jtr1w                                                   
      do  730  ji=1,mxsi                                                
      jiw=(ji+li-2)*2-isiw                                              
      if(jiw.lt.jimin) go to 730                                        
      if(jiw.gt.jimax) go to 730                                        
      call ovlapa(ovlpi,lf,jf,li,ji,ic,ktzf3)                           
      if(ktout(2).eq.0) go to 540                                       
      if(iskip1) 522,522,524                                            
  522 iskip1=1                                                          
      call skip(iskip2,npgs,numrun,ititol)                              
  524 l1 =liw/2                                                         
      fj1=float(jiw)/2.0                                                
      l3 =lfw/2                                                         
      fj3=float(jfw)/2.0                                                
      iskip2=iskip2+1                                                   
c     *********************************************                     
c     if(l1-l3.ne.4) ovlpi=(0.d0,0.d0)                                  
c     *********************************************                     
      if(iskip2.gt.noline) call skip(iskip2,npgs,numrun,ititol)         
*     write(6,530) l1,fj1,l3,fj3,ovlpi                                  
c 530 format(1h ,2x,8hone step,6x,3hli=,i3,3x,3hji=,f5.1,23x,3hlf=,i3,  
c    13x,3hjf=,f5.1,5x,7hovlap=(,2e15.7,2h ))                           
  530 format(1h ,2x,i3,3x,f5.1,3x,i3,3x,f5.1,5x,2e15.7)                 
  540 fac7=fac3*u9(lfw,isfw,ltr1w,is1w,jfw,jtr1w,liw,isiw,jiw)          
     1         *ovlpi                                                   
      if(ktzf3.gt.0) fac7=fac7*frfac(1,nch)                             
      do 600 msi1=1,mxsi                                                
      msiw=2*msi1-isiw-2                                                
      fac8=fac7*cleb(liw,0,isiw,msiw,jiw,msiw)                          
      do 600 msf1=1,mxsf                                                
      msfw=2*msf1-isfw-2                                                
      do 590 mjr1=1,mxjp                                                
      mlfw=2*mjr1+mshif+msfw-msiw                                       
      if(mlfw)550,560,560                                               
  550 mlfwa=-mlfw                                                       
      phas=(-1.0)**(mlfwa/2)                                            
      go to 565                                                         
  560 mlfwa=mlfw                                                        
      phas=1.0                                                          
  565 mlgd=mlfwa/2+1                                                    
      if(mlgd-mmax) 570,570,590                                         
  570 n=(lfw-mlfwa)/2+1                                                 
      k=(lfw+mlfwa)/2+1                                                 
      mlfwp=-mlfw                                                       
      mjfwp=-mlfw+msfw                                                  
      mjrwp=msiw+mlfw-msfw                                              
      fac9=fac8*phas*exp((faclog(n)-faclog(k))/2.0)                     
     1         *cleb(lfw,mlfwp,isfw ,msfw ,jfw,mjfwp)                   
     2         *cleb(jfw,mjfwp,jtr1w,mjrwp,jiw,msiw )                   
      do 580 j=1,jmax                                                   
      betn=fac9*(plmd(j,mlgd)+(0.d0,2.d0)/pie*plmq(j,mlgd))*beta(2)     
      betf=fac9*(plmd(j,mlgd)-(0.d0,2.d0)/pie*plmq(j,mlgd))*beta(2)     
  580 samp(j,mjr1,msf1,msi1)=samp(j,mjr1,msf1,msi1)+(betn+betf)/2.d0    
*     samp(j,mjr1,msf1,msi1)=samp(j,mjr1,msf1,msi1)+betn/2.d0           
* 580 samf(j,mjr1,msf1,msi1)=samf(j,mjr1,msf1,msi1)+betf/2.d0           
  590 continue                                                          
  600 continue                                                          
  730 continue                                                          
  710 continue                                                          
  800 continue                                                          
  810 continue                                                          
  820 continue                                                          
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine writea                                                 
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/pot/vd(2),wd(2),vsod(2),wsod(2),rrd(2),ard(2),rcd(2),      
     1csdgd(2),rid(2),aid(2),rsord(2),asord(2),rsoid(2),asoid(2),       
     2ktop(2),rgd(2),agd(2),ktisp(2),risrd(2),aisrd(2),visd(2),risid(2),
     3aisid(2),wisd(2),csd,drhod(12),rhobn,rhoin,rhosor,rhosoi,rhoisr,  
     4rhoisi,rhobc,rhobng,fkayar,fkayai,fkasor,fkasoi,fkaisr,fkaisi,    
     5fkaybg,vbe,wbe,vsoe,wsoe,vise,wise,vbcd(905,2),wbcd(905,2),       
     6vsov(905,2),wsov(905,2),pnloc(2)                                  
      common/coul/z(16),sigld(90,2),fd(90,2),gd(90,2),fpd(90,2),        
     1gpd(90,2)                                                         
      call skip(k,npgs,numrun,ititol)                                   
*     write(6,3000)                                                     
 3000 format(1h ,5x,46h*****(   basic data for every channel   )*****/) 
      write(6,3050) elabi                                               
 3050 format(1h ,5x,16hincident channel,8x,7h elabi=,f9.5)              
      i=1                                                               
   20 write(6,3100) ecmd(i),qvlue(i),fk(i),z(i),drd(i),pnloc(i),        
     1              pmas(i),pspn(i),pz(i),tmas(i),tspn(i),tz(i),        
     2              vd(i),wd(i),vsod(i),wsod(i),visd(i),wisd(i),        
     3              rrd(i),rid(i),rsord(i),rsoid(i),risrd(i),risid(i),  
     4              ard(i),aid(i),asord(i),asoid(i),aisrd(i),aisid(i),  
     5              rcd(i),rgd(i),agd(i),csdgd(i),lmind(i),lmaxd(i)     
 3100 format( 7h   ecm=,f11.5,2x,7h qvleu=,                             
     1 f11.5,2x,7h     k=,f9.5,4x,7h   eta=,f9.5,4x,                    
     1 7h    dr=,f9.5,4x,7h pnloc=,f9.5,/,                              
     1 7h  pmas=,f9.5,4x,7h  pspn=,f9.5,4x,                             
     2 7h    pz=,f9.5,4x,7h  tmas=,f9.5,4x,                             
     2 7h  tspn=,f9.5,4x,7h    tz=,f9.5,/,                              
     3 7h     v=,f9.5,4x,7h     w=,f9.5,4x,                             
     3 7h   vso=,f9.5,4x,7h   wso=,f9.5,4x,                             
     4 7h   vis=,f9.5,4x,7h   wis=,f9.5,/,                              
     4 7h    r0=,f9.5,4x,7h    ri=,f9.5,4x,                             
     5 7h  rsor=,f9.5,4x,7h  rsoi=,f9.5,4x,                             
     5 7h  risr=,f9.5,4x,7h  risi=,f9.5,/,                              
     6 7h    ar=,f9.5,4x,7h    ai=,f9.5,4x,                             
     6 7h  asor=,f9.5,4x,7h  asoi=,f9.5,4x,                             
     7 7h  aisr=,f9.5,4x,7h  aisi=,f9.5,/,                              
     7 7h    rc=,f9.5,4x,7h    rg=,f9.5,4x,                             
     8 7h    ag=,f9.5,4x,7h  csdg=,f9.5,4x,                             
     8 7h  lmin=,i5  ,8x,7h  lmax=,i5)                                
      go to (30,60),i                                                   
   30 write(6,3150)                                                     
 3150 format(1h ,5x,12hexit channel)                                  
      i=2                                                               
      go to 20                                                          
   60 write(6,3250) rmax,drd(1),nrmin,nrmax                             
 3250 format(17h integration data,5x,7h  rmax=,f9.5,4x,                 
     17h    dr=,f9.5,4x,7h nrmin=,i5,8x,7h nrmax=,i5,/)                
*     write(6,3260) ttrj                                                
 3260 format(22h total transfered spin,4x,7h jtrns=,f6.2)                
*     write(6,3300) beta(2)                                             
 3300 format(//22h channel mixing factor,8x,8hone step,5x,f16.5)        
      return                                                            
      end                                                               
c     ******************************************************************
      subroutine writeb                                                 
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/ffc/amp,d02,fts,ta,tb,fnrng,ktff(5),nrange,nafac,nbfac,lsa 
     1,lsb,iread(15),jjj,fread(28),wread(21),ktrl(20),ff(905),ffi(905), 
     2ffone(905),nrangd(2,4),n1xd(2,4),itap(66,10,2),frfac(2,4),frso(2) 
      complex*16 ffone,frfac                                            
      common/zrff/fis,fls,fjs,ffout(10),bb,bsubl,nrn,lss,redwd,zeta,    
     1fma,fmas,fmb,fmbs,sa,fja,sb,fjb,za,zas,zb,zbs,vb(905)             
      integer idata(11)                                                 
      data idata/4h  ga,4huss ,4h yuk,4hawa ,4hinit,4hial ,4hfina       
     1                 ,4hl   ,4hj1  ,4h    ,4hj2  /                    
    4 write(6,3330) fjs,lss,fis,fts                                     
 3330 format(17h transfered spins,3x,6hjspin=,f5.2,3x,6hlspin=,i3,3x,     
     1  6hsspin=,f5.2,3x,6htspin=,f5.2)                                  
      if(ktff(1).eq.8) go to 5                                          
*     write(6,3333) amp,ktff(4),d02,fnrng                               
 3333 format(1h ,5x,22hmixing amplitude  amp=,f9.5,15x,8hktff(4)=,i3,   
     19x,4hd02=,-4pf11.5,4he+04,11x,6hfnrng=,0pf9.5,/)                  
*   5 write(6,3335)ktff(1),ktff(2),ktff(3),dr                           
    5 continue
 3335 format(23h form factor parameters,//,12h    ktff(1)=,i3,4x,       
     19h ktff(2)=,i3,4x,9h ktff(3)=,i3,10x,3hdr=,f9.5)                  
   20 kt=ktff(1)                                                        
      go to (25,30,35,40,90,45,100,110),kt                              
   25 continue
*  25 write(6,3340)(fread(i),i=1,5),bsubl,(fread(i),i=6,11)             
 3340 format(1h ,9x,2hv=,f9.5,4x,3hr0=,f9.5,4x,3har=,f9.5,4x,3hrc=,f9.5,
     +3x,4hrce=,
     1f9.5,7h bsubl=,f9.2/1h ,9x,2hw=,f9.5,4x,3hri=,f9.5,4x,3hai=,
     +f9.5,4x,3hrg=,f9.5,5x,2hb=,f9.5,4x,3hcs=,f9.5/)   
      go to 50      
   30 write(6,3342)(fread(i),i=1,4),iread(1),lss,fread(5),nrn,zeta,redwd
     1,bb,fread(6),fread(7)                                             
 3342 format(1h ,9x,2he=,f9.5,4x,3hzz=,f9.5,4x,3hm1=,f9.5,4x,3hm2=,
     +f9.5,5x,2hn=,i3,
     111x,2hl=,i3/1h ,8x,3hrn=,f9.5,3x,4hnrn=,i6,2x,8hff(nrn)=,f9.5,
     +7h redwd=,f9.5,2x,5hbeta=,f9.5,2x,5h-b/a=,f9.5,5x,2ha=,f9.5/)  
      go to 50                                                          
   35 write(6,3344)fread(1),iread(1),lss                                
 3344 format(1h ,8x,3hnu=,f9.5,5x,2hn=,i3,11x,2hl=,i3/)                 
      go to 50                                                          
*  40 write(6,3345) fread(23)                                           
   40 continue
 3345 format(1h ,85x,6hpnloc=,f9.5)                                      
      write(6,3346)fread(1),iread(2),fread(2),fread(3),iread(1),lss,fjs,
     1ffout(6),fread(6),fread(8),fread(7),fread(9),frso(1),frso(2)
 3346 format(1h ,3x,2he=,f9.5,4x,3hzz=,i3,10x,3hm1=,f9.5,4x,3hm2=,f9.5,
     +3x,4hnod=,i3,9x,2hl=,i3,11x,2hj=,f5.2/1h ,
     +3x,2hv=,f9.5,4x,3hr0=,f9.5,4x,3har=,f9.5,4x,3hrc=,f9.5/1h ,
     +3x,4hvso=,f9.5,2x,4hrso=,f9.5,3x,4haso=,f9.5/)                    
      go to 50                                                          
   90 write(6,3347)                                                     
 3347 format(1h ,2x,17hform factor input/)                              
      go to 50                                                          
   45 if(iread(3))46,46,47                                              
   46 itype1=idata(1)                                                   
      itype2=idata(2)                                                   
      go to 48                                                          
   47 itype1=idata(3)                                                   
      itype2=idata(4)                                                   
   48 itype3=idata(5)                                                   
      itype4=idata(6)                                                   
      itype5=idata(7)                                                   
      itype6=idata(8)                                                   
      write(6,3348)itype1,itype2,fread(10),fread(11),fis,fts            
 3348 format(1h ,10x,2a4,16htype interaction,6x,3hv0=,f9.5,
     +7h range=,f9.5,5x,
     12hs=,f5.2,9x,2ht=,f5.2)                                           
   49 write(6,3350)fread(2),fread(3),itype3,itype4,fread(1),iread(9),   
     1fread(21),fread(9),iread(1),iread(4),fread(12),fread(6),fread(8), 
     2fread(7),itype5,itype6,fread(13),iread(10),fread(22),fread(20),   
     3iread(5),iread(6),fread(14),fread(17),fread(19),fread(18)         
 3350 format(1h ,8x,3hm1=,f9.5,4x,3hm2=,f9.5/1h ,2a4,1x,2he=,f9.5,4x,
     +3hzz=,i3,11x,
     12hv=,f9.5,3x,4hvso=,f9.5,3x,4hnod=,i3,11x,2hl=,i3,11x,2hj=,
     +f5.2/1h ,8x,3hr0=,
     2f9.5,4x,3har=,f9.5,4x,3hrc=,f9.5/1h ,2a4,1x,2he=,f9.5,4x,3hzz=,
     +i3,11x,2hv=,
     3f9.5,3x,4hvso=,f9.5,3x,4hnod=,i3,11x,2hl=,i3,11x,2hj=,f5.2/1h ,
     +8x,3hr0=,f9.5
     4,4x,3har=,f9.5,4x,3hrc=,f9.5/)                                    
      go to 50                                                          
  100 itype3=idata(9)                                                   
      itype4=idata(10)                                                  
      itype5=idata(11)                                                  
      itype6=idata(10)                                                  
      write(6,3352)fread(10),fread(11),fis,fts                          
 3352 format(1h ,9x,2hk=,f9.5,3x,4heta=,f9.5,5x,2hs=,f5.2,9x,2ht=,
     +f5.2)         
      go to 49                                                          
  110 write(6,3354) (fread(i),i=1,12)                                   
 3354 format(1h ,6x,5hpmas=,f9.5,2x,5htmas=,f9.5,4x,3hpz=,f4.0,9x,
     +3htz=,f4.0,   
     17x,5hpspn=,f6.2,5x,5htspn=,f6.2/1h ,8x,3hs1=,f6.2,7x,3hl1=,
     +f6.2,7x,3hj1=,f6.2,
     27x,3hs2=,f6.2,7x,3hl2=,f6.2,7x,3hj2=,f6.2/)                       
   50 if(ktff(5).eq.0.or.ktff(5).eq.2) go to 80                         
   75 call skip(k,npgs,numrun,ititol)                                   
      write(6,3360) nrmax                                               
 3360 format(1h ,37hform factor used for overlap integral,3x,5hnr=1,i3/)
      write(6,3362)                                                     
 3362 format(1h ,9hreal part/)                                          
      write           (6,3365)(nr,ff(nr),nr=1,nrmax)                    
 3365 format(5(1h ,i6,e16.8))                                           
      if((ktff(1).ne.1).and.(ktff(1).ne.5).and.(ktff(1).ne.8)) go to 80 
      if(ktff(3).eq.0) go to 80                                         
      if(nrmax.gt.130) call skip(k,npgs,numrun,ititol)                  
      write(6,3370)                                                     
 3370 format(1h ,15h imaginary part/)                                   
      write(6,3365)(nr,ffi(nr),nr=1,nrmax)                              
   80 if(ktff(5).lt.2) go to 5300                                       
      fnc0=10.0                                                         
      fktjls=ktff(4)                                                    
      write(7,3375) fnc0,fts,amp,fktjls                                 
 3375 format(8f10.5)                                                    
      fnc=10.53                                                         
      write(7,3380) fnc                                                 
 3380 format(f7.2/8h(5e16.8))                                           
      write(7,3390) (ff(nr),nr=1,nrmax)                                 
 3390 format(5e16.8)                                                    
      if(ktff(3).eq.0) go to 5300                                       
      fnci=10.54                                                        
      write(7,3380) fnci                                                
      write(7,3390) (ffi(nr),nr=1,nrmax)                                
 5300 continue                                                          
 6000 continue                                                          
 7000 return                                                            
      end                                                               
c     ******************************************************************
      subroutine xsm                                                    
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      common/theta/jmax,j1min,j1max,ktheta,mmax,nnnn,
     1thetad(181),sint(181),    
     1cost(181),plmd(181,10),xdwb(181),poldw(181),asymdw(181),
     2tenpol(181,6),plmq(181,10)                      
      common/bamp/samp(181,10,3,3),samf(181,10,3,3) 
      common/xsecs/xdwl                            
      complex*16 samp,tb,samf,rmval,ucn,thamp,samx(181,10,3,3)         
      dimension a(3),b(3,3),xdwl(15,181)
      character*12 ddsn   
      common/dataset/ddsn                                               
      data a/0.,2.,1.414213562/,b/0.7071067814,-1.414213562,0.7071067814
     1      ,1.224744871,-1.224744871,0.,1.732050808,0.,0./ 
      ucn=(0.d0,1.d0)
      do 30 j=1,jmax                                                    
       xdwb(j)=0.0 
       do k=1,15  
        xdwl(k,j)=0.0  
       enddo                                                  
       do 20 k=1,6                                                      
        tenpol(j,k)=0.0                                                 
   20 continue                                                          
      do 30 msi1=1,mxsi                                                 
      do 30 msf1=1,mxsf                                                 
      do 30 mjr1=1,mxjp                                                 
      mj=2*mjr1+mshif                                                   
      if(mj) 30,22,23                                                   
   22 fac=1.0                                                           
      go to 24                                                          
   23 fac=2.0                                                           
   24 xdwb(j)=xdwb(j)+fac*samp(j,mjr1,msf1,msi1)                        
     1             *conjg(samp(j,mjr1,msf1,msi1))                       
   30 continue                                                          
      fac=0.06332574*fk(2)*(2.0*tspn(2)+1.0)                            
     1   /(ecmd(1)*ecmd(2)*fk(1)*(2.0*pspn(1)+1.0)*(2.0*tspn(1)+1.0))   
      do 35 j=1,jmax                                                    
   35 xdwb(j)=fac*xdwb(j)                                               
      if(pspn(2).lt.0.1) go to 180                                      
      do 100 j=1,jmax                                                   
      t0=0.0                                                            
      t1=0.0                                                            
      t20=0.                                                            
      t21=0.                                                            
      t22=0.                                                            
      do 95 ma=1,mxsi                                                   
      do 95 mb=1,mxsf                                                   
      do 95 mjp=1,mxjp                                                  
      mj=2*mjp+mshif                                                    
      if(mj) 90,82,83                                                   
   82 fac1=1.0                                                          
      go to 85                                                          
   83 fac1=2.0                                                          
   85 ta=fac1*samp(j,mjp,mb,ma)*conjg(samp(j,mjp,mb,ma))                
      t0=t0+ta                                                          
      t20=t20+b(mb,1)*ta                                                
      if(mb.eq.mxsf) go to 90                                           
      tb=fac1*samp(j,mjp,mb,ma)*conjg(samp(j,mjp,mb+1,ma))              
      t1=t1+a(mxsf)*dimag(tb)                                           
      t21=t21+b(mb,2)*tb                                                
      if(mb.gt.1) go to 90                                              
      t22=t22+fac1*b(mb,3)*samp(j,mjp,mb,ma)*conjg(samp(j,mjp,mb+2,ma)) 
   90 continue                                                          
   95 continue                                                          
      if(t0.lt.1.d-60) t0=1.d+35                                        
      poldw(j)=100.0*t1/t0                                              
      if(mxsf.lt.3) go to 100                                           
      tenpol(j,1)=100.*t22/t0                                           
      tenpol(j,2)=100.*t21/t0                                           
      tenpol(j,3)=100.*t20/t0                                           
  100 continue                                                          
  180 if(pspn(1).lt.0.1) go to 205                                      
      do 200 j=1,jmax                                                   
      t0=0.0                                                            
      t1=0.0                                                            
      t20=0.                                                            
      t21=0.                                                            
      t22=0.                                                            
      do 195 ma=1,mxsi                                                  
      do 195 mb=1,mxsf                                                  
      do 195 mjp=1,mxjp                                                 
      mj=2*mjp+mshif                                                    
      if(mj) 190,182,183                                                
  182 fac1=1.0                                                          
      go to 185                                                         
  183 fac1=2.0                                                          
  185 ta=fac1*samp(j,mjp,mb,ma)*conjg(samp(j,mjp,mb,ma))                
      t0=t0+ta                                                          
      t20=t20+b(ma,1)*ta                                                
      if(ma.eq.mxsi) go to 190                                          
      tb=fac1*samp(j,mjp,mb,ma)*conjg(samp(j,mjp,mb,ma+1))              
      t1=t1+a(mxsi)*dimag(tb)                                           
      t21=t21+b(ma,2)*tb                                                
      if(ma.gt.1) go to 190                                             
      t22=t22+fac1*b(ma,3)*samp(j,mjp,mb,ma)*conjg(samp(j,mjp,mb,ma+2)) 
  190 continue                                                          
  195 continue                                                          
      if(t0.lt.1.d-60) t0=1.d+35                                        
      asymdw(j)=-50.0*t1/t0*b(1,3)                                      
      if(mxsi.lt.3) go to 200                                           
      tenpol(j,4)=100.*t22/t0                                           
      tenpol(j,5)=100.*t21/t0                                           
      tenpol(j,6)=100.*t20/t0                                           
  200 continue                                                          
  205 continue 
*     -------------------------------------------------------------------  
      if(ktout(1).ne.8) goto 7000
*     -------------------------------------------------------------------
*     read Euler angles in case of amplitudes being rotated (ktout(1)=8)
*     -------------------------------------------------------------------   
      if(mxjp.gt.5) then
       print*,' Too many m-values for samp definitions -so redimension'
       stop
      endif  
*     ------------------------------------------------------------------- 
*     read the Euler angles from front6z generated file
      open(49,file='rot.'//ddsn,status='unknown') 
      read(49,*) alpha,betaa,gamaa
      alpha=3.1415926536d0*alpha
      betaa=3.1415926536d0*betaa
      gamaa=3.1415926536d0*gamaa 
      print*,' Euler angles have been read from file ' 
      print*
      print*,real(alpha),real(betaa),real(gamaa)
      print*
*     ------------------------------------------------------------------- 
*     copy the original m>0 matrix elements upward into samf array  
*     anticipating using 2*mxjp elements - the lowest mxjp for m<0
      do j=1,jmax
       do msi1=1,mxsi                                      
        do msf1=1,mxsf                                                  
         do mjr1=1,mxjp 
         samf(j,mjr1+mxjp,msf1,msi1)=samp(j,mjr1,msf1,msi1)       
         enddo
        enddo
       enddo
      enddo    
*     ------------------------------------------------------------------- 
*     retrieve the l and j transfer and spin values 
      rjval=trj(1,2) 
      lval=ltr(1,2)
      rsi=pspn(1)
      rsf=pspn(2)
*     print*,rsi,rsf
*     print*,trs(1,2),ltr(1,2),trj(1,2) 
*     ------------------------------------------------------------------- 
*     compute all m<0 elements using the Madison symmetry relation       
      do j=1,jmax
       do msi1=1,mxsi   
        rmsi=msi1-rsi-1.d0  
        msip=nint(-rmsi+rsi+1)                                    
        do msf1=1,mxsf   
         rmsf=msf1-rsf-1.d0 
         msfp=nint(-rmsf+rsf+1)                                         
         do mjr1=1,mxjp 
          rmf1=0.5d0*(2*mjr1+mshif)
          mjrp=mxjp+1-mjr1  
          phase=(-1)**nint(lval+rsi-rmsi+rsf-rmsf+rjval+rmf1)
          samf(j,mjrp,msf1,msi1)=phase*samp(j,mjr1,msfp,msip)       
         enddo
        enddo
       enddo
      enddo
*     this completes the calculation of the negative m amplitudes  
*     samf stores the full t-matrix in the Madison coordimates
*     ------------------------------------------------------------------- 
*     zero samp, prior to rotating amplitudes samf (into samp) 
*     ------------------------------------------------------------------- 
      do msi1=1,mxsi                                       
       do msf1=1,mxsf                                                   
        do mjr1=1,2*mxjp 
         do j=1,jmax          
          samx(j,mjr1,msf1,msi1)=(0.d0,0.d0)     
         enddo
        enddo
       enddo
      enddo    
*     ------------------------------------------------------------------- 
*     rotate the samf for the (j,m) to other coordinate system
*     ------------------------------------------------------------------- 
      do msi1=1,mxsi                                       
       do msf1=1,mxsf                                                   
        do mjr1=1,2*mxjp
         rmf1=mjr1-rjval-1.d0 
         do mjrp=1,2*mxjp
          rmfp=mjrp-rjval-1.d0  
*         compute the star of the rotation matrix
          rmval=rmat(rjval,rmfp,rmf1,betaa)
          rmval=rmval*exp(ucn*(rmfp*alpha+rmf1*gamaa))
          do j=1,jmax                   
           samx(j,mjr1,msf1,msi1)=samx(j,mjr1,msf1,msi1)
     *                     +rmval*samf(j,mjrp,msf1,msi1) 
          enddo      
         enddo
        enddo
       enddo
      enddo
*     ------------------------------------------------------------------- 
*     diagnostic print May 2012
*     ------------------------------------------------------------------- 
*     do msi1=1,mxsi                                       
*     do msf1=1,mxsf                                                   
*     do mjr1=1,2*mxjp
*     if(abs(samx(10,mjr1,msf1,msi1)).gt.0.1d-6) then
*     print*,samx(10,mjr1,msf1,msi1)-samf(10,mjr1,msf1,msi1) 
*     endif
*     enddo
*     enddo
*     enddo
*     ------------------------------------------------------------------- 
*     cross sections, and also the sum for checking purposes
*     ------------------------------------------------------------------- 
      do mjr1=1,2*mxjp          
       do j=1,jmax
        do msi1=1,mxsi                                                 
         do msf1=1,mxsf                                                 
          thamp=samx(j,mjr1,msf1,msi1)                                 
          xdwl(mjr1,j)=xdwl(mjr1,j)+thamp*conjg(thamp)                 
         enddo  
        enddo   
*       compute cross section summed over m-substates         
        xdwl(15,j)=xdwl(15,j)+xdwl(mjr1,j)        
       enddo 
      enddo                                                            
*     ------------------------------------------------------------------- 
*     multiply by phase space factor 
*     ------------------------------------------------------------------- 
      do mjr1=1,15
       do j=1,jmax                                                    
        xdwl(mjr1,j)=fac*xdwl(mjr1,j)       
       enddo
      enddo 
      open(20,file='20.'//ddsn,status='unknown') 
      open(19,file='19.'//ddsn,status='unknown')
*     ----------------------------------------------------------------
      write(20,3199)(thetad(j),xdwl(15,j),j=1,jmax)    
 3199 format(0pf7.2,1pe15.6)   
      write(19,3198)(thetad(j),(xdwl(k,j),k=1,9),xdwl(15,j),j=1,jmax)   
 3198 format(f7.2,9e11.4,e12.5)
*     -------------------------------------------------------------------                                                                 
 7000 return                                                            
      end                                                               
c     ******************************************************************
      subroutine plm(xsect,ndata,jmax,thetad)                           
      implicit real*8(a-h,o-z)                                          
      common/base/ijt,nubchn,beta(2),nrmin,nrmax,nr3min,nr3max,rmax,dr, 
     1drd(2),elabi,qvlue(2),ecm,ecmd(2),fkay,fk(2),eta,drho,pspn(2),    
     2pz(2),pmas(2),tspn(2),tz(2),tmas(2),zz(2),tautau(2),fmud(2),lmax, 
     3ll,lmind(2),lmaxd(2),trs(2,2),ltr(2,2),trj(2,2),ic,lc,fjc,ichnl,  
     4kcheck,npgs,numrun(5),noline,ititol(12),iskip1,iskip2,ktout(10),  
     5isiw,isfw,jpw,mxsi,mxsf,mxjp,mmmm,ttrj,mshif                      
      dimension xsect(181),flet(101),ix(5),thetad(181),fno(10)          
      data blank,star/1h ,1h*/                                          
      data fno/1h ,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1h /                 
      smax=xsect(1)                                                     
      do 10 i=1,ndata                                                   
      do 10 j=1,jmax                                                    
      if(xsect(j).gt.smax) smax=xsect(j)                                
   10 continue                                                          
      smax=log10(smax)                                                  
      imax=int(smax)                                                    
      if(smax.gt.0.0) imax=imax+1                                       
      do 20 i=1,5                                                       
   20 ix(i)=imax-5+i                                                    
      fix=ix(1)                                                         
      call skip(k,npgs,numrun,ititol)                                   
*     write(6,100) ix                                                   
  100 format(1h ,//,i15,4i25/11x,2h10,23x,2h10,23x,2h10,23x,2h10,23x,
     +2h10/    
     113x,101h*------------------------*------------------------*-------
     2-----------------*------------------------*)                      
      do 70 j=1,jmax                                                    
      do 30 i=1,101                                                     
   30 flet(i)=blank                                                     
      do 40 i=1,ndata                                                   
      ii=(log10(xsect(j))-fix)*25.0+1.0                                 
      if(ndata.eq.1) go to 50                                           
      sq=fno(i)                                                         
      if(ii.lt.1) go to 40                                              
      if(flet(ii).ne.blank) sq=star                                     
      flet(ii)=sq                                                       
   40 continue                                                          
      go to 60                                                          
   50 if(ii.lt.1) go to 60                                              
      flet(ii)=star                                                     
   60 continue
*  60 write(6,200) thetad(j),flet                                       
  200 format(f11.2,2x,101a1)                                            
   70 continue                                                          
*     write(6,300)                                                      
  300 format(13x,101h*------------------------*------------------------*
     1------------------------*------------------------*)               
      return                                                            
      end                                                               
*************************************************************************
      subroutine sim(fa,res,m,n,h)
*------------------------------------------------------------------------
*     subroutine does the integral of fa stored
*     in the arrays of the same name using simpsons rule. the step length
*     is h and the integral is between the elements m and n of the arrays
*     only. resulting integral is placed in res.
*------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension fa(*),dq(905)
      do 90 i=m,n
      dq(i)=fa(i)
   90 continue
      rq1=dq(m+1)
      rq2=dq(m+2)
      i=m+3
   98 continue
      if(i.ge.n) go to 99
      rq1=rq1+dq(i)
      rq2=rq2+dq(i+1)
      i=i+2
      go to 98
   99 continue
      res=0.33333333333d0*h*(dq(m)+4.d0*rq1+2.d0*rq2-dq(n))
      return
      end
c     ******************************************************************
      real*8 function rmat(aj,am,an,beta)
      implicit real*8(a-h,o-z)
c     ------------------------------------------------------------------
c     small d rotation matrices for non-integer j. brink and satchler
c     formula. rotation of beta radians about y-axis.
c     ------------------------------------------------------------------
      common/clebma/faclog(500)
      integer t
      t=0
      rmat=1.d0
      if(abs(am).gt.aj.or.abs(an).gt.aj) go to 888
      if(aj.eq.0) return
      sc=sin(beta/2.d0)
      dc=cos(beta/2.d0)
      ind1=nint(aj+am+1.d0)
      ind2=nint(aj-an+1.d0)
      ind3=nint(an-am+1.d0)
      ind4=nint(aj-am+1.d0)
      ind5=nint(aj+an+1.d0)
      rnum=(faclog(ind1)+faclog(ind4)+faclog(ind2)+faclog(ind5))
      rnum=rnum/2.d0
      rmat=0.d0
      go to 997
   88 ind6=nint(2.d0*(aj-t)+am-an)
      ind7=nint(2.d0*t+an-am)
      den=faclog(ind1)+faclog(ind2)+faclog(ind3)+faclog(t+1)
      rmat=rmat+(-1.d0)**t*dexp(rnum-den)*dc**(ind6)
     1     *sc**(ind7)
   99 t=t+1
      ind1=ind1-1
      ind2=ind2-1
      ind3=ind3+1
  997 if((ind1.le.0).or.(ind2.le.0)) return
      if(ind3.le.0) go to 99
      go to 88
  888 rmat=0.d0
      return
      end
      SUBROUTINE COULFG(XX,ETA1,XLMIN,XLMAX, FC,GC,FCP,GCP,
     *                  MODE1,KFN,IFAIL)
C  REVISED #5 IJT WITH L-T ALGORITHMN FOR CONTINUED FRACTIONS,
C  AND IFAIL > 0 FOR AVOIDED EXPONENT CHECKS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           C
C                                                                      C
C  A. R. BARNETT           MANCHESTER  MARCH   1981                    C
C                                                                      C
C  ORIGINAL PROGRAM 'RCWFN'      IN    CPC  8 (1974) 377-395           C
C                 + 'RCWFF'      IN    CPC 11 (1976) 141-142           C
C  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314           C
C  THIS VERSION WRITTEN UP       IN    CPC 27 (1982) 147-166           C
C                                                                      C
C  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA1 (INCLUDING 0), C
C   AND REAL LAMBDA(XLMIN) .GT. -1 FOR INTEGER-SPACED LAMBDA VALUES    C
C   THUS GIVING POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER   C
C   EQUATION,TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF     C
C   THE DIRAC EQUATION ,ALSO SPHERICAL & CYLINDRICAL BESSEL EQUATIONS  C
C                                                                      C
C  FOR A RANGE OF LAMBDA VALUES (XLMAX - XLMIN) MUST BE AN INTEGER,    C
C  STARTING ARRAY ELEMENT IS M1 = MAX0(IDINT(XLMIN+ACCUR),0) + 1       C
C      SEE TEXT FOR MODIFICATIONS FOR INTEGER L-VALUES                 C
C                                                                      C
C  IF 'MODE' = 1  GET F,G,F',G'   FOR INTEGER-SPACED LAMBDA VALUES     C
C            = 2      F,G      UNUSED ARRAYS MUST BE DIMENSIONED IN    C
C            = 3      F               CALL TO AT LEAST LENGTH (1)      C
C  IF 'KFN'  = 0 REAL        COULOMB FUNCTIONS ARE RETURNED            C
C            = 1 SPHERICAL   BESSEL      "      "     "                C
C            = 2 CYLINDRICAL BESSEL      "      "     "                C
C  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT                          C
C                                                                      C
C  PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF 'MACHINE ACCURACY'    C
C   IN OSCILLATING REGION X .GE. ETA1 + SQRT(ETA1**2 + XLM(XLM+1))     C
C   COULFG IS CODED FOR REAL*8 ON IBM OR EQUIVALENT  ACCUR = 10**-16   C
C   USE AUTODBL + EXTENDED PRECISION ON HX COMPILER  ACCUR = 10**-33   C
C   FOR MANTISSAS OF 56 & 112 BITS. FOR SINGLE PRECISION CDC (48 BITS) C
C   REASSIGN DSQRT=SQRT ETC.  SEE TEXT FOR COMPLEX ARITHMETIC VERSION  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION    FC(1),GC(1),FCP(1),GCP(1)
      LOGICAL      ETANE0,XLTURN
      COMMON       /STEE / PACCQ,NFP,NPQ,IEXP,M1
C***  COMMON BLOCK IS FOR INFORMATION ONLY.  NOT REQUIRED IN CODE
C***  COULFG HAS CALLS TO: DSQRT,DABS,DMOD,IDINT,DSIGN,DFLOAT,DMIN1
      DATA ZERO,ONE,TWO,TEN2,ABORT /0.0D0, 1.0D0, 2.0D0, 1.0D2, 2.0D4/
      DATA HALF,TM30,BIG / 0.5D0, 1.0D-30, 1.0D+100/
      DATA RT2DPI /0.79788 45608 02865 35587 98921 19868 76373 D0/
C *** THIS CONSTANT IS  DSQRT(TWO/PI):  USE Q0 FOR IBM REAL*16: D0 FOR
C ***  REAL*8 & CDC DOUBLE P:  E0 FOR CDC SINGLE P; AND TRUNCATE VALUE.
C
                        ACCUR = 1.0D-16
C ***            CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
      MODE  = 1
      IF(MODE1 .EQ. 2 .OR. MODE1 .EQ. 3 ) MODE = MODE1
      IFAIL = 0
      IEXP  = 1
      NPQ   = 0
      ETA   = ETA1
      GJWKB = ZERO
      PACCQ = ONE
      IF(KFN .NE. 0) ETA = ZERO
                 ETANE0  = ETA .NE. ZERO
      ACC   = ACCUR * 10D0
      ACC4  = ACC*TEN2*TEN2
      ACCH  = DSQRT(ACC)
C ***    TEST RANGE OF XX, EXIT IF.LE.DSQRT(ACCUR) OR IF NEGATIVE
C
      IF(XX .LE. ACCH)                          GO TO 100
      X     = XX
      XLM   = XLMIN
      IF(KFN .EQ. 2)  XLM = XLM - HALF
      IF(XLM .LE. -ONE .OR. XLMAX .LT. XLMIN)   GO TO 105
      E2MM1 = ETA*ETA + XLM*XLM + XLM
      XLTURN= X*(X - TWO*ETA) .LT. XLM*XLM + XLM
      DELL  = XLMAX - XLMIN + ACC
      IF(DABS(DMOD(DELL,ONE)) .GT. ACC) WRITE(6,2040)XLMAX,XLMIN,DELL
      LXTRA = IDINT(DELL)
      XLL   = XLM + DFLOAT(LXTRA)
C ***       LXTRA IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
C ***       XLL  IS MAX LAMBDA VALUE, OR 0.5 SMALLER FOR J,Y BESSELS
C ***         DETERMINE STARTING ARRAY ELEMENT (M1) FROM XLMIN
      M1  = MAX0(IDINT(XLMIN + ACC),0) + 1
      L1  = M1 + LXTRA
C
C ***    EVALUATE CF1  =  F   =  FPRIME(XL,ETA,X)/F(XL,ETA,X)
C
      XI  = ONE/X
      FCL = ONE
      PK  = XLL + ONE
      PX  = PK  + ABORT
      F   =  ETA/PK + PK*XI
         IF(DABS(F).LT.TM30) F = TM30
         D = ZERO
         C = F
C
C ***   BEGIN CF1 LOOP ON PK = K = LAMBDA + 1
C
    4 PK1   = PK + ONE
        EK  = ETA / PK
        RK2 = ONE + EK*EK
        TK  = (PK + PK1)*(XI + EK/PK1)
        D   =  TK - RK2 * D
        C   =  TK - RK2 / C
         IF(DABS(C).LT.TM30) C = TM30
         IF(DABS(D).LT.TM30) D = TM30
         D = ONE/D
         DF = D * C
         F  = F * DF
            IF(D .LT. ZERO) FCL = - FCL
         PK = PK1
                          IF( PK .GT. PX ) GO TO 110
      IF(DABS(DF-ONE) .GE. ACC)             GO TO 4
                  NFP = PK - XLL - 1
      IF(LXTRA .EQ. 0)                          GO TO 7
C
C *** DOWNWARD RECURRENCE TO LAMBDA = XLM. ARRAY GC,IF PRESENT,STORES RL
C
      FCL = FCL/BIG
      FPL = FCL*F
      IF(MODE .EQ. 1) FCP(L1) = FPL
                      FC (L1) = FCL
      XL  = XLL
      RL  = ONE
      EL  = ZERO
      DO 6  LP = 1,LXTRA
         IF(ETANE0) EL = ETA/XL
         IF(ETANE0) RL = DSQRT(ONE + EL*EL)
         SL    =  EL  + XL*XI
         L     =  L1  - LP
         FCL1  = (FCL *SL + FPL)/RL
         FPL   =  FCL1*SL - FCL *RL
         FCL   =  FCL1
         FC(L) =  FCL
         IF(MODE .EQ. 1) FCP(L)  = FPL
         IF(MODE .NE. 3 .AND. ETANE0) GC(L+1) = RL
         if(abs(FCL).gt.BIG) then
                do 55 LP1=L,M1+LXTRA
                IF(MODE .EQ. 1) FCP(LP1) = FCP(LP1)*1d-20
 55                                     FC (LP1) = FC(LP1)*1d-20
                FCL=FC(L)
                FPL=FPL*1d-20
                endif
    6 XL = XL - ONE
      IF(FCL .EQ. ZERO) FCL = ACC
      F  = FPL/FCL
C ***    NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
C ***    EVALUATE CF2 = P + I.Q  AGAIN USING STEED'S ALGORITHM
C ***    SEE TEXT FOR COMPACT COMPLEX CODE FOR SP CDC OR NON-ANSI IBM
C
    7 IF( XLTURN ) CALL JWKB(X,ETA,DMAX1(XLM,ZERO),FJWKB,GJWKB,IEXP)
      IF( IEXP .GT. 1 .OR. GJWKB .GT. ONE/(ACCH*TEN2))  GO TO 9
          XLTURN = .FALSE.
      TA =  TWO*ABORT
      PK =  ZERO
      WI =  ETA + ETA
      P  =  ZERO
      Q  =  ONE - ETA*XI
      AR = -E2MM1
      AI =  ETA
      BR =  TWO*(X - ETA)
      BI =  TWO
      DR =  BR/(BR*BR + BI*BI)
      DI = -BI/(BR*BR + BI*BI)
      DP = -XI*(AR*DI + AI*DR)
      DQ =  XI*(AR*DR - AI*DI)
    8 P     = P  + DP
         Q  = Q  + DQ
         PK = PK + TWO
         AR = AR + PK
         AI = AI + WI
         BI = BI + TWO
         D  = AR*DR - AI*DI + BR
         DI = AI*DR + AR*DI + BI
         C  = ONE/(D*D + DI*DI)
         DR =  C*D
         DI = -C*DI
         A  = BR*DR - BI*DI - ONE
         B  = BI*DR + BR*DI
         C  = DP*A  - DQ*B
         DQ = DP*B  + DQ*A
         DP = C
         IF(PK .GT. TA)                         GO TO 120
      IF(DABS(DP)+DABS(DQ).GE.(DABS(P)+DABS(Q))*ACC)   GO TO 8
                      NPQ   = PK/TWO
                      PACCQ = HALF*ACC/DMIN1(DABS(Q),ONE)
                      IF(DABS(P) .GT. DABS(Q)) PACCQ = PACCQ*DABS(P)
C
C *** SOLVE FOR FCM = F AT LAMBDA = XLM,THEN FIND NORM FACTOR W=W/FCM
C
      GAM = (F - P)/Q
            IF(Q .LE. ACC4*DABS(P))             GO TO 130
      W   = ONE/DSQRT((F - P)*GAM + Q)
            GO TO 10
C *** ARRIVE HERE IF G(XLM) .GT. 10**6 OR IEXP .GT. 70 & XLTURN = .TRUE.
    9 W   = FJWKB
      GAM = GJWKB*W
      P   = F
      Q   = ONE
C
C *** NORMALISE FOR SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
C
   10                     ALPHA = ZERO
          IF(KFN  .EQ. 1) ALPHA = XI
          IF(KFN  .EQ. 2) ALPHA = XI*HALF
                          BETA  = ONE
          IF(KFN  .EQ. 1) BETA  = XI
          IF(KFN  .EQ. 2) BETA  = DSQRT(XI)*RT2DPI
      FCM  = DSIGN(W,FCL)*BETA
           FC(M1)  = FCM
                      IF(MODE .EQ. 3)           GO TO 11
           IF(.NOT. XLTURN)   GCL =  FCM*GAM
           IF(      XLTURN)   GCL =  GJWKB*BETA
           IF( KFN .NE. 0 )   GCL = -GCL
           GC(M1)  = GCL
           GPL =  GCL*(P - Q/GAM) - ALPHA*GCL
                      IF(MODE .EQ. 2)           GO TO 11
           GCP(M1) = GPL
           FCP(M1) = FCM*(F - ALPHA)
   11 IF(LXTRA .EQ. 0 ) RETURN
C *** UPWARD RECURRENCE FROM GC(M1),GCP(M1)  STORED VALUE IS RL
C *** RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
C ***    XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
         W    = BETA*W/DABS(FCL)
         MAXL = L1 - 1
      DO 12 L = M1,MAXL
                      IF(MODE .EQ. 3)           GO TO 12
                      XL = XL + ONE
         IF(ETANE0)   EL = ETA/XL
         IF(ETANE0)   RL = GC(L+1)
                      SL = EL + XL*XI
         GCL1     = ((SL - ALPHA)*GCL - GPL)/RL
                IF(ABS(GCL1).gt.BIG) GO TO 140
         GPL      =   RL*GCL -  (SL + ALPHA)*GCL1
         GCL      = GCL1
         GC(L+1)  = GCL1
                      IF(MODE .EQ. 2)           GO TO 12
         GCP(L+1) = GPL
         FCP(L+1) = W*(FCP(L+1) - ALPHA*FC(L+1))
   12 FC(L+1)     = W* FC(L+1)
      RETURN
C
C ***    ERROR MESSAGES
C
  100 IFAIL = -1
      WRITE(6,2000) XX,ACCH
 2000 FORMAT(' FOR XX = ',1P,D12.3,' TRY SMALL-X  SOLUTIONS',
     *' OR X NEGATIVE'/ ,' SQUARE ROOT ACCURACY PARAMETER =  ',D12.3/)
      RETURN
  105 IFAIL = -2
      WRITE (6,2005) XLMAX,XLMIN,XLM
 2005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES:XLMAX,XLMIN,XLM = ',
     *1P,3D15.6/)
      RETURN
  110 IFAIL = -3
      WRITE (6,2010) ABORT,F ,DF,PK,PX,ACC
 2010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',F10.0,' ITERATIONS',/
     *' F,DF,PK,PX,ACCUR =  ',1P,5D12.3//)
      RETURN
  120 IFAIL = -4
      WRITE (6,2020) ABORT,P,Q,DP,DQ,ACC
 2020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',F7.0,' ITERATIONS',/
     *' P,Q,DP,DQ,ACCUR =  ',1P,4D17.7,D12.3//)
      RETURN
  130 IFAIL = -5
      WRITE (6,2030) P,Q,ACC,DELL,LXTRA,M1
 2030 FORMAT(' FINAL Q.LE.DABS(P)*ACC*10**4 , P,Q,ACC = ',1P,3D12.3,4X,
     *' DELL,LXTRA,M1 = ',D12.3,2I5 /)
      RETURN
 2040 FORMAT(' XLMAX - XLMIN = DELL NOT AN INTEGER ',1P,3D20.10/)
  140 IFAIL=L1-L
      RETURN
      END
C
      SUBROUTINE JWKB(XX,ETA1,XL,FJWKB,GJWKB,IEXP)
      REAL*8          XX,ETA1,XL,FJWKB,GJWKB,DZERO
C *** COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS    FOR XL.GE. 0
C *** AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
C *** CALLS DMAX1,SQRT,ALOG,EXP,ATAN2,FLOAT,INT        BARNETT FEB 1981
      DATA   ZERO,HALF,ONE,SIX,TEN/ 0.0E0, 0.5E0, 1.0E0, 6.0E0, 10.0E0 /
      DATA  DZERO, RL35, ALOGE  /0.0D0, 35.0E0, 0.43429 45 E0 /
      X     = XX
      ETA   = ETA1
      GH2   = X*(ETA + ETA - X)
      XLL1  = DMAX1(XL*XL + XL,DZERO)
      IF(GH2 + XLL1 .LE. ZERO) RETURN
       HLL  = XLL1 + SIX/RL35
       HL   = SQRT(HLL)
       SL   = ETA/HL + HL/X
       RL2  = ONE + ETA*ETA/HLL
       GH   = SQRT(GH2 + HLL)/X
       PHI  = X*GH - HALF*( HL*ALOG((GH + SL)**2/RL2) - ALOG(GH) )
          IF(ETA .NE. ZERO) PHI = PHI - ETA*ATAN2(X*GH,X - ETA)
      PHI10 = -PHI*ALOGE
      IEXP  =  INT(PHI10)
      IF(IEXP .GT. 70) GJWKB = TEN**(PHI10 - FLOAT(IEXP))
      IF(IEXP .LE. 70) GJWKB = EXP(-PHI)
      IF(IEXP .LE. 70) IEXP  = 0
      FJWKB = HALF/(GH*GJWKB)
      RETURN
      END
