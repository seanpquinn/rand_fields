c 15  Jan 2003: fixed whicharm bug (c.f. TWJL email of 12 Dec '02)
c 24 June 2002: added calculations of path lengths through LISM components
c 26  May 2002: modified for NE2001 routines which are cleaned-up
c             versions of development routines.
c Nov 1999 - May 2002: development versions
c 1992-1993: TC93 version

c--CHANGE-substituted smtau (not used by us anyways)
c--with sungcdistance. This is an input comming from
c--hammurabi.
c-------------------------------------------------------
      subroutine dmdsm(
     .   l,b,ndir,dmpsr,dist,limit,sm,sungcdistance,smtheta,smiso)
      implicit none
      real l,b,dmpsr,dist,sm,sungcdistance,smtheta,smiso
c-------------------------------------------------------
      integer ndir
      character*1 limit

c  Computes pulsar distance and scattering measure
c  from model of Galactic electron distribution.

c  Input: real l        galactic longitude in radians
c         real b        galactic latitude in radians
c         integer ndir  >= 0 calculates dist from dmpsr
c                       < 0 for dmpsr from dist
c Input or output:
c         real dmpsr    (dispersion measure in pc/cm^3)
c         real dist     (distance in kpc)

c  Output:
c         char*1 limit  (set to '>' if only a lower distance limit can be
c                        given; otherwise set to ' ')
c         sm            (scattering measure, uniform weighting) (kpc/m^{20/3})
c         smtau         (scattering measure, weighting for pulse broadening)
c         smtheta       (scattering measure, weighting for angular broadening
c                        of galactic sources)
c         smiso         (scattering measure appropriate for calculating the
c                       isoplanatic angle at the source's location
c       parameter(alpha = 11./3.)
c       parameter(pi = 3.14159)
c       parameter(c_sm = (alpha - 3.) / 2. * (2.*pi)**(4.-alpha) )


        real c_sm, c_u, sm_factor
        parameter(c_sm = 0.181)         ! constant in sm definition
        parameter(c_u = 10.16)          ! units conversion for sm
        parameter(sm_factor = c_sm * c_u)


        integer wg1, wg2, wga, wggc, wglism, wgcN, wgvN
        common /modelflags/ wg1, wg2, wga, wggc, wglism, wgcN, wgvN

c---------------------------------------------------------------------------------------------
c flag for printing ne at given position to fort.26

        integer wfat
c---------------------------------------------------------------------------------------------



c parameters of large-scale components (inner+outer+arm components):
        real n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
        common/galparams/n1h1,h1,A1,F1,n2,h2,A2,F2,
     .                na,ha,wa,Aa,Fa

c factors for controlling individual spiral arms:
c       narm:   multiplies electron density (in addition to the`fac'
c                     quantities)
c       warm:   arm width factors that multiply nominal arm width
c       harm:   arm scale height factors
c       farm:   factors that multiply n_e^2 when calculating SM

        integer narmsmax, narmsmax1
        parameter (narmsmax=5, narmsmax1=narmsmax+1)
        real narm, warm, harm, farm
        common/armfactors/
     .     harm(narmsmax),narm(narmsmax),warm(narmsmax),farm(narmsmax)

        real armpaths, armdistances
        common/armpathlengths/ armpaths(narmsmax1), armdistances(narmsmax1)




        real dx0, dy0, dz0
        common/dxyz/dx0,dy0,dz0

        integer whicharm

c Large scale components:

        real ne1, ne2, nea
        real F1val, F2val, Faval

c Galactic center:

        real negc, Fgc

c LISM:
        real nelism, Flism
        integer wlism, wLDR, wLHB, wLSB, wLOOPI
        real ldr_path, lhb_path, lsb_path, loopI_path
        real ldr_dist, lhb_dist, lsb_dist, loopI_dist
        integer wtemp

c clumps:
        real necN, FcN
        integer hitclump

c voids:
        real nevN, FvN
        integer hitvoid, wvoid

c subroutines needed:
c       density_2001 (and those that it calls) in density.NE2001.f
c       scattering routines in scattering98.f

        real R0, rrmax, zmax, dmax
c------------------------------------------------------------------
c--CHANGE  Substitued R0 by sungcdistance
c       data R0/8.5/
c------------------------------------------------------------------
c       data rrmax/30.0/                ! Max radius for reliable ne
        data rrmax/50.0/                ! Max radius for reliable ne
c       data zmax/1.76/                 ! Max |z|
c       data zmax/5.00/                 ! Max |z|
        data zmax/25.00/                ! Max |z|
        data dmax/50.0/                 ! maximum distance calculated

        logical first
        data first/.true./

        save

c other variables
        real x, y, z, r, rr
        real sl, cl, sb, cb

        real d, dstep, dtest, dstep_pc, dd

        real dm, dmstep
        real sm_sum1, sm_sum2, sm_sum3, sm_sum4, sm_term
        real sm_sum1_last, sm_sum2_last, sm_sum3_last, sm_sum4_last
        integer nstep
c       integer ncount
        integer i

        real dm1, dm2, dma, dmgc, dmlism, dmcN, dmvN
        real sm1, sm2, sma, smgc, smlism, smcN, smvN
        real dsm1, dsm2, dsma, dsmgc, dsmlism, dsmcN, dsmvN

        integer wtotal
        real ne

c------------------------------------------------------------------------
c do not think the if below alters l,b,dist

c--CHANGE--------------------------------------------------------------------
        R0=sungcdistance
c------------------------------------------------------------------



        if(first) then
c initial call to density routine to set variable values
c through read-in of parameter file:
        x = 0.0
        y = R0
        z = 0.0
            call density_2001(x,y,z,
     .        ne1,ne2,nea,negc,nelism,necN,nevN,
     .        F1val, F2val, Faval, Fgc, Flism, FcN, FvN,
     .        whicharm, wlism, wldr, wlhb, wlsb, wloopI,
     .        hitclump, hitvoid, wvoid, sungcdistance)
        first=.false.
        endif
c-------------------------------------------------------------------------
        sl=sin(l)
        cl=cos(l)
        sb=sin(b)
        cb=cos(b)
        limit=' '
c       dstep=0.02                      ! Step size in kpc
c       dstep = min(h1, h2) / 10.       ! step size in terms of scale heights
        dstep=0.01
        if(ndir.lt.0) dtest=dist
        if(ndir.ge.0) dtest=dmpsr/(n1h1/h1)   ! approximate test distance
        nstep = dtest / dstep           ! approximate number of steps
        if(nstep.lt.10) dstep=dtest/10  ! make # steps >= 10

c  Sum until dm is reached (ndir >= 0) or dist is reached (ndir < 0).
c  Guard against too few terms by counting number of terms (ncount) so that
c  routine will work for n_e models with large n_e near the Sun.

    5   continue
        dstep_pc = 1000.*dstep
        dm=0.0
        sm_sum1 = 0.                    ! sum of C_n^2
        sm_sum2 = 0.                    ! sum of C_n^2 * s
        sm_sum3 = 0.                    ! sum of C_n^2 * s^2
        sm_sum4 = 0.                    ! sum of C_n^2 * s^{5./3.}

        do i=1,narmsmax1
          armpaths(i) = 0.
          armdistances(i) = 0.
        enddo

        dm1 = 0.
        dm2 = 0.
        dma = 0.
        dmgc = 0.
        dmlism = 0.
        dmcN = 0.
        dmvN = 0.

        sm1 = 0.
        sm2 = 0.
        sma = 0.
        smgc = 0.
        smlism = 0.
        smcN = 0.
        smvN = 0.

        ldr_path = 0.
        lhb_path = 0.
        lsb_path = 0.
        loopI_path = 0.

        ldr_dist = 0.
        lhb_dist = 0.
        lsb_dist = 0.
        loopI_dist = 0.

c        ncount = 0

        continue



        d=-0.5*dstep

c--------------------------CHANGE-28feb2005--------------------------------------------------AA
c       do 10 i=1,99999
c          ncount = ncount + 1
          d=dist                        ! Distance from Sun in kpc
c--------------------------------------------------------------------------------------------------
c exchanging d by dist

          r=dist*cb
          z=dist*sb

c         r=d*cb
          x=r*sl
          y=R0-r*cl
c         z=d*sb

c--------------------------------------------------------------------------------------------------


          rr=sqrt(x**2 + y**2)          ! Galactocentric radius

c---------------------------------------------------------------------------------------------------
c call to density, why? --> I dont get the lt. 3, but it is necessary to call density eventualy

          if(ndir.lt.3) then
            call density_2001(x,y,z,
     .        ne1,ne2,nea,negc,nelism,necN,nevN,
     .        F1val, F2val, Faval, Fgc, Flism, FcN, FvN,
     .        whicharm, wlism, wldr, wlhb, wlsb, wloopI,
     .        hitclump, hitvoid, wvoid,sungcdistance)
          endif
c---------------------------------------------------------------------------------------------------

          ne=
     .       (1.-wglism*wlism)*
     .       (wg1*ne1 +
     .        wg2*ne2 +
     .        wga*nea +
     .        wggc*negc) +
     .        wglism*wlism*nelism
          ne = (1-wgvN*wvoid)*ne + wgvN*wvoid*nevN + wgcN*necN

                                sm=ne

        return
        end
