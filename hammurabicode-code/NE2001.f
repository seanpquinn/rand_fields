c 	program NE2001

c--------------------------------------------------------------------------
c
c should obtain the ne density in the file fort.32. If ne without LISM is desired, change destination filename, and see density.NE2001.f
c where some lines should be erased to avoid calling the LISM routines
c
c--------------------------------------------------------------------------


c  calls dmdsm to get pulsar distance as a function of l, b, and DM,
c  or pulsar DM as a function of l, b, and distance.

	implicit none
	integer nargs
	character*80 inbuf
c---------------------------------------------------------------------------
c originaly 3 real and one integer value. now they are corresponding arrays
	integer ndata
	integer j
	parameter (ndata=500000)
	real aldeg(ndata), abdeg(ndata), admd(ndata)
	integer andir(ndata)
	real ldeg, bdeg, dmd
	integer ndir
c---------------------------------------------------------------------------
	real dist
	real dm, sm, smtau, smtheta, smiso
	real tau, sbw, stime, theta_g, theta_x, transfreq, emsm 

	real nu
        data nu/1./

	character*1 limit

	real rad
	data rad/57.2957795/

	real vperp
        data vperp/100./

c functions:
c-----------------------------------------------------------------------------
c altered zone below

c	integer iargc
c	external iargc
        real em, tauiss, scintbw, scintime, theta_xgal, theta_gal
	real transition_frequency

c	nargs = iargc()
c	if(nargs .ge. 1) then

c	for generating the magfield ne template (i.e. ne without the LISM), the file to which we save is noLISM.32, instead of fort.32
	open(34, file='fort.32', status='unknown')
 	open(33, file='cppinput.inp', status='old')
	do j=1,ndata
	   read (33,*) aldeg(j), abdeg(j), admd(j), andir(j)
	   ldeg=aldeg(j)
	      bdeg=abdeg(j)
c dmd is the distance from earth to source in kpc
	      dmd=admd(j)
	      ndir=-1
c	   call getarg(1, inbuf)
c	   read(inbuf, *) ldeg
c	   call getarg(2, inbuf)
c	   read(inbuf, *) bdeg
c	   call getarg(3, inbuf)
c	   read(inbuf, *) dmd
c	   call getarg(4, inbuf)
c          read(inbuf, *) ndir
	
c           close (11)
c----------------------------------------------------------------------

c	else
c	   write(*,*) 'Usage: NE2001 l b DM/D ndir'
c           write(*,*) '       l (deg)'
c           write(*,*) '       b (deg)'
c           write(*,*) '       DM/D (pc cm^{-3} or kpc)'
c           write(*,*) '       ndir = 1 (DM->D)   -1 (D->DM) '
c	   stop
c        endif

c	write(*,*)'#NE2001 input: 4 parameters'
c	write(*,"(f10.4, t20, a, t30, a, t55, a)") 
c     .        ldeg, 'l', '(deg)', 'GalacticLongitude'
c	write(*,"(f10.4, t20, a, t30, a, t55, a)") 
c     .        bdeg, 'b', '(deg)', 'GalacticLatitude'
c	write(*,"(f10.4, t20, a, t30, a, t55, a)") 
c     .        dmd,  'DM/D', '(pc-cm^{-3}_or_kpc)', 
c     .        'Input_DM_or_Distance'
c	write(*,"(i10, t20, a, t30, a, t55, a)") 
c     .        ndir, 'ndir', '1:DM->D;-1:D->DM', 'Which?(DM_or_D)'
c	write(*,*)'#NE2001 output: 14 values'
c
	if(ndir.ge.0) then
          dm = dmd
c-----------------------------------------------------------------------------
c notice that dmdsm is called with ldeg/rad, and bdeg/rad
c therefore l,b values in dmdsm are all in rads.
c-----------------------------------------------------------------------------
	  call dmdsm(ldeg/rad,bdeg/rad,ndir, dm,
     .               dist,limit,sm,smtau,smtheta,
     .               smiso)

c	  write(*,"(a,f9.4, t20, a, t30, a, t55,a)") 
c     .         limit, dist, 'DIST', '(kpc)','ModelDistance'
c	  write(*,"(f10.4, t20, a, t30, a, t55, a)") 
c     .         dm, 'DM', '(pc-cm^{-3})', 'DispersionMeasure' 
c	  write(*,"(f10.4, t20, a, t30, a, t55, a)") 
c     .         dm*abs(sin(bdeg/rad)), 'DMz', 
c     .            '(pc-cm^{-3})', 'DM_Zcomponent' 
	else
	  dist=dmd
	  call dmdsm(ldeg/rad,bdeg/rad,ndir, dm,
     .               dist,limit,sm,smtau,smtheta,
     .               smiso)
c	  write(6,*) 'dm,sm,smtau,smtheta = ', 
c    .               dm,sm,smtau,smtheta
c1020	  format(f8.2,3(1x,e8.3))
c	  write(6,*) 'dmz = ', dm*abs(sin(bdeg/rad))
c	  write(*,"(f10.4, t20, a, t30, a, t55, a)") 
c     .          dist, 'DIST', '(kpc)', 'Distance'
c	  write(*,"(f10.4, t20, a, t30, a, t55,a)") dm, 'DM', 
c     .            '(pc-cm^{-3})','ModelDM' 
c	  write(*,"(f10.4, t20, a, t30, a, t55, a)") 
c     .            dm*abs(sin(bdeg/rad)), 'DMz', 
c     .            '(pc-cm^{-3})', 'model' 
	endif
c-----------------------------------------------------------------------------------------
c the calculation of scattering parameters is unecessary for our purposes
	end do
	close (34)
	close (33)
c-----------------------------------------------------------------------------------------
c calculate scattering parameters

 
        tau = tauiss(dist, smtau, nu)
        sbw = scintbw(dist, smtau, nu)
        stime = scintime(smtau, nu, vperp)
        theta_x = theta_xgal(sm, nu)
        theta_g = theta_gal(smtheta, nu)
	transfreq = transition_frequency(sm,smtau,smtheta,dist)
        emsm = em(sm)
	stop
	end
