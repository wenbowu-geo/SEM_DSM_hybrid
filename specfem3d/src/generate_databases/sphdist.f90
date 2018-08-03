!*-----        Location Routines      ------------------- + ----
!                                                       ydist
      subroutine distazbaz(dlats,dlons,dlatr,dlonr,delta,cazim,bazim)
!
!     AUTHOR:  Brian L.N. Kennett  RSES, ANU
!     DATE:    January 1985
!     PURPOSE:
!             YDIST        Calculates distance and azimuth
!                          for spheroidal earth between
!                          specified geographic source and
!                          receiver station coordinates
!
!     USAGE:ydist dlons dlats dlonr dlonr
!     OUTPUT: gcarc backazimuth azimuth
!*---------------------------------------------------*------*
!     PARAMETERS
!
      implicit none
      double precision dlats, dlons, dlatr, dlonr, delta, cazim, bazim
!
!     dlats  latitude of source
!     dlons  longitude of source
!     dlatr  latitude of receiver
!     dlonr  longitude of receiver
!     delta  angular distance
!     cazim  apparent azimuth at an array
!     bazim   azimuth from epicentre to receiver
!
!------------------------------------------------------*------*
!
!     implicit real*8 (a-h,o-z)
      double precision  ecc,re,ec1,pi,pib2,degr,rlats,rlons,rlatr
      double precision  rlonr,glats,glatr,sps,cps,spr,cpr,rs,rr
      double precision  trs,prs,trr,prr,AS,BS,CS,DS,ES,GS,HS,KS
      double precision  AR,BR,CR,DR,ER,GR,HR,KR
      double precision  cosdr,deltar,sindr,deltak,szs,czs,szr,czr
      double precision  azima
!      double precision  gra
!      character(len=40) temp
!                          radius on spheroid
!      gra(x,y,e) = dsqrt( (1.0d0-e)**2 /
!     &                   ((1.0d0-e*y)**2 + e*e*x*y ) )
!      !ecc = 0.000001
!      !re = 6378.388
      ecc=0.0
      re=6371
      ec1 = (1.0d0-ecc)**2
      pi = 3.141592653589793
      pib2 = pi/2.0
      degr = pi/180.0
      
!      CALL getArg(1,temp)
!      read(temp,*)dlons
!      CALL getArg(2,temp)
!      read(temp,*)dlats
!      CALL getArg(3,temp)
!      read(temp,*)dlonr
!      CALL getArg(4,temp)
!      read(temp,*)dlatr

      rlats = dlats*degr
      rlons = dlons*degr
      rlatr = dlatr*degr
      rlonr = dlonr*degr
!                          geocentric coordinates
      glats = datan2 ( ec1*dsin(rlats) ,dcos(rlats) )
      glatr = datan2 ( ec1*dsin(rlatr) ,dcos(rlatr) )
      sps = dsin(glats)**2
      cps = dcos(glats)**2
      spr = dsin(glatr)**2
      cpr = dcos(glatr)**2
!                          radii at source,receiver
      rs = re*gra(sps,cps,ecc)
      rr = re*gra(spr,cpr,ecc)
!
      trs = pib2 - glats
      prs = dlons*degr
      trr = pib2 - glatr
      prr = dlonr*degr
!                          direction cosines for source
      AS = dsin(trs)*dcos(prs)
      BS = dsin(trs)*dsin(prs)
      CS = dcos(trs)
      DS = dsin(prs)
      ES = -dcos(prs)
      GS = dcos(trs)*dcos(prs)
      HS = dcos(trs)*dsin(prs)
      KS = -dsin(trs)
!                          direction cosines for receiver
      AR = dsin(trr)*dcos(prr)
      BR = dsin(trr)*dsin(prr)
      CR = dcos(trr)
      DR = dsin(prr)
      ER = -dcos(prr)
      GR = dcos(trr)*dcos(prr)
      HR = dcos(trr)*dsin(prr)
      KR = -dsin(trr)
!                          distance
      cosdr = AS*AR + BS*BR + CS*CR
      deltar = dacos(cosdr)
      sindr = dsin(deltar)
!
      deltak = deltar*0.5d0*(rr+rs)
      delta = deltar/degr
!                          azimuth
      szs = DS*AR + ES*BR
      czs = GS*AR + HS*BR + KS*CR
      szr = DR*AS + ER*BS
      czr = GR*AS + HR*BS + KR*CS
!                          azima - azimuth to source
!                          bazim - backazimuth from source
!                          cazim - apparent azimuth at an array
      if (szr.eq.0.0) then
        bazim = 0.0
        if(dlats.gt.dlatr)then
           azima = 360.0
        else
           azima = 180.0
        endif
      else
        bazim = datan2(-szs ,-czs ) /degr
        azima = datan2(-szr ,-czr ) /degr
      end if
      if( bazim .lt. 0.0) bazim = bazim + 360.0
      cazim = azima + 180.0
      if( azima .lt. 0.0) azima = azima + 360.0
!
      if( cazim.lt. 0.0) cazim = cazim + 360.0
!
!      write(*,*)delta, bazim, azima
      contains
      function gra(x,y,e)
        double precision x,y,e
        double precision gra
          gra = dsqrt( (1.0d0-e)**2 / ((1.0d0-e*y)**2 + e*e*x*y ) )
      end function gra
      end subroutine
