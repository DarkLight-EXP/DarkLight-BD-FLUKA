*$ CREATE MAGFLD.FOR
*COPY MAGFLD
*
*===magfld=============================================================*
*
      !SUBROUTINE MAGFLD ( X, Y, Z, BTX, BTY, BTZ, B)
      SUBROUTINE MAGFLD ( X, Y, Z, T, BTX, BTY, BTZ, B, NREG, IDISC)

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'

*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1988-2009      by Alberto Fasso` & Alfredo Ferrari *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     Created  in 1988    by     Alberto Fasso`, CERN - TIS            *
*                                                                      *
*     Last change on 15-oct-09     by    Alfredo Ferrari               *
*                                                                      *
*     Input variables:                                                 *
*            x,y,z = current position                                  *
*            nreg  = current region                                    *
*     Output variables:                                                *
*            btx,bty,btz = cosines of the magn. field vector           *
*            B = magnetic field intensity (Tesla)                      *
*            idisc = set to 1 if the particle has to be discarded      *
*                                                                      *
*----------------------------------------------------------------------*
*
*     Last change on 12-sep-12    by    Advanced FLUKA course teacher  *
*
      ! INCLUDE '(RTDFCM)'
      ! INCLUDE '(LTCLCM)'

*     lfirst, to flag first call of routine
*      LOGICAL LFIRST
*      DATA    LFIRST / .TRUE. /
*      SAVE    LFIRST

*     magnetic data:
*     - dipole field (T):
      ! DATA DIPFLD / 8.34D+00 /
      ! SAVE DIPFLD

*     magnetic regions:
      ! SAVE NUMQUA, NUMDIP

*     default values
      ! IDISC = 0
      ! BTX = ZERZER
      ! BTY = ZERZER
      ! BTZ = ONEONE
      ! B   = ZERZER

      BTX = 0.
      BTY = 0.
      BTZ = 0.
      B   = 0.

*----------------------------------------------------------------------*
*     first call: initialisation
*----------------------------------------------------------------------*
      !DL target sits at (-465, 76, 2941.5)
      
      STRK1 = 0.62 ! strength of first permanent magnet - 0.084*QM20*length/aperture QM20=2, aper=2.6, len = 9.3
      STRK2 = 1. ! strength of second permanent magnet - 0.084*QM21*length/aperture QM21=2.8
      STRK3 = 0.62 ! strength of third permanent magnet - 0.084*QM22*length/aperture QM22=1.7
      STRK4 = 0.0 ! strenght of original permanent magnet quad
      STRK5 = 0.22 ! 0.012705*QM18*length/aperture current in beam optics model of 8.636, aper = 3.3, len = 10
      STRK6 = 0.12    ! 0.012705*QM19*length/aperture current in beam optics model of -6.986, aper = 3.3, len = 10
      !STRK1 = -0.382 !0.012705*QM5*length/aperture current in beam optics model of -9.915 A, aper = 3.3, len = 10
      !STRK2 = 0.258  !0.012705*QM6*length/aperture current in beam optics model of 6.704  A, aper = 3.3, len = 10

      APR = 2.6 !cm for PMQs
      APR2 = 3.3 !cm for EMQs  !!
 
      x_local = x + 465. !subtract extra 20cm for the amount the target has shifted downstream
      y_local = y - 76.
      z_local = z - 2941.5

      x1 = 20. + 9.3/2 ! centre of PMQ1
      x2 = x1 + 9.3/2 + 3. + 9.3/2.! centre of PMQ2
      x3 = x2 + 9.3/2 + 3.5 + 9.3/2.! centre of PMQ3
      x4 = x3 + 9.3/2 + 20.  + 9.3/2 ! centre of PMQ4
      x5 = x4 + 9.3/2 + 3.25 + 9.3/2.! centre of PMQ5
      x6 = x5 + 9.3/2 + 5. + 3.211 + 0.539 + 7.5 + 5.! centre of EMQ5
      x7 = x6 + 5. + 3.5 + 5.!centre of EMQ6

      !print*,x1,x2,x3,x4,x5,x6,x7

      !new fake distances to check FLUKA without overlapping fields for the same strenghts - target out
      !x1 = 20. + 9.3/2
      !x2 = x1 + 9.3/2 + 13. + 9.3/2
      !x3 = x2 + 9.3/2 + 13. + 9.3/2
      !x4 = x3 + 9.3/2 + 13. + 9.3/2
      !x5 = x4 + 9.3/2 + 13. + 9.3/2
      !x6 = x5 + 9.3/2 + 10.911 + 5.
      !x7 = x6 +5. +10. +5.

      !Aveen here is when you gonna do

      CALL PMQUAD(-STRK1,APR,x1,x_local,y_local,z_local,BTX,BTY,BTZ)
      CALL PMQUAD(STRK2,APR,x2,x_local,y_local,z_local,BTX,BTY,BTZ)
      CALL PMQUAD(-STRK3,APR,x3,x_local,y_local,z_local,BTX,BTY,BTZ)
      CALL PMQUAD(STRK4,APR,x4,x_local,y_local,z_local,BTX,BTY,BTZ)
      CALL PMQUAD(-STRK4,APR,x5,x_local,y_local,z_local,BTX,BTY,BTZ)
      CALL MSQ(STRK5,APR2,x6,x_local,y_local,z_local,BTX,BTY,BTZ)
      CALL MSQ(-STRK6,APR2,x7,x_local,y_local,z_local,BTX,BTY,BTZ)

*     get field module:
      B   = SQRT( BTX**2 + BTY**2 + BTZ**2 )

*     check against numerical precision
      IF ( B .GT. 1.0D-12 ) THEN
*        normalise field components - should be direction cosines
*     cosine of the angle of the vector dotted with coordinate over magnitude
         BTX = BTX / B
         BTY = BTY / B
         BTZ = BTZ / B
      !WRITE(713,*) x_local,y_local,z_local,BTX*B,BTY*B,BTZ*B,B
      ENDIF
      !B=B/1.e6 
      RETURN
*=== End of subroutine magfld =========================================*
      END

      SUBROUTINE MSQUAD(STRK,APR,Z0,x_,y_,z_,Bx,By,Bz)
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'

      ! INCLUDE '(DBLPRC)'
      ! INCLUDE '(DIMPAR)'
      ! INCLUDE '(IOUNIT)'

c     Increments Bx,By,Bz
c     3D field distribution from a short quadrupole, based on Baartman R. Quadrupole shapes. Physical Review Special Topics-Accelerators and Beams. 2012 Jul 30;15(7):074002.  
c     INPUTS:
c     STRK Integrated field gradient in TESLA
c     APR  aperture ~ pole tipe radius IN CM
c     Z0   position of the quad center along the beam axis
c     x_,y_,z_ coordinates at which the magnetic field is evaluted (in cm)
c     INCREMENTED:
c     Bx,By,Bz the 3 component of the magnetic field in TESLA
      PI = 4.D0*DATAN(1.D0)
      al = 4./PI*APR !length scale for the sech^2 k(z) on-axis gradient function
      !rotate the local coordinate system to account for the fact that magnetic quads are rotated +45 deg wrt to electro static quads (see PSAB paper referred to above)
      x = 1./SQRT(2.0)*( x_ - y_)
      y = 1./SQRT(2.0)*( x_ + y_)
      ! choosing del x and del y at same time then we might get Bz

      !if beam is outside of the aperture, 
      !or too far away from the magnet (>6*aperture, where the field should have dropped by nearly 4 orders of magnitue, good enough), 
      !return zero magnetic field
      if((x**2+y**2.GT.APR**2).OR.(ABS(z_-z0).GE.(6*APR))) then
            Bx = 0.0
            By = 0.0
            Bz = 0.0
            return
      endif

      !normalize x,y,z by the characteristic length l, and center the field around z0
      x = x/al
      y = y/al
      z = (z_-z0)/al

      !now the field components, copy-pasted directly from mathematica
      Bx_loc = (-0.5*STRK*(Sin(2*x))/(Cos(2*x) + Cosh(2*z)))
      By_loc = (0.5*STRK*(Sin(2*y))/(Cos(2*y) + Cosh(2*z)))
      
      !now rotate back Bx and By -45 deg
      Bx = Bx + 1./SQRT(2.0)*( Bx_loc + By_loc)
      By = By + 1./SQRT(2.0)*(-Bx_loc + By_loc)

      !Increment Bz as is, since the rotation is done around z 
      Bz = Bz - 0.5*STRK*Sinh(2*z)*(1/(Cos(2*y)+Cosh(2*z))
     *- 1/(Cos(2*x)+Cosh(2*z)))
      
      return
      end



      SUBROUTINE MSQ(STRK,APR,X0,x_,y_,z_,Bx,By,Bz)
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'

      ! INCLUDE '(DBLPRC)'
      ! INCLUDE '(DIMPAR)'
      ! INCLUDE '(IOUNIT)'
c     Wrapper subroutine to change from Baartman's coordinate system to Fluka
      
      !initialize variables
            !why are these initialized to zero?
      Bx_l = 0.
      By_l = 0.
      Bz_l = 0.

      !Call previous subroutine with swapped coordinates
      CALL MSQUAD(STRK,APR,-X0,z_,y_,-x_,Bz_l,By_l,Bx_l)

      !increment new coordinates
      Bx = Bx - Bx_l
      Bz = Bz + Bz_l
      By = By + By_l

      !print*,Bx,By,Bz
      !initial Bx = Bx - Bz_l and Bz = Bz + Bx_l does not work, onyly changes sign of Bx
      !changed to Bx = Bz - Bz_l and Bz = Bx + Bx_l which gives z the correct field with wrong sign
      !but Bx still incorrect, should return zero but its returning same as Bz...why?
      !dont think I can call Bx and Bz within one another, its making them both the same
      !trying Bx = Bx - Bx_l and Bz = Bz + Bz_l seems to work better? swaps z and x succesfully!
      return
      End

      SUBROUTINE PMQUAD(STRK,APR,X0,x_,y_,z_,Bx,By,Bz)
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
      !wrapper around MSQUAD to model the DarkLight permament magnet quads. Kinda kludgy, but it will do
c     INPUTS:
c     X0   position of the quad center along the beam axis
c     x_,y_,z_ coordinates at which the magnetic field is evaluted (in cm)
c     OUTPUTS:
c     Bx,By,Bz the 3 component of the magnetic field in TESLA
      !APR=5.2/2.0 !in cm
      PI=4.D0*DATAN(1.D0)
      al = 4/PI*APR
      alen = al*0.658479 !good enough approximation for al*ACosh(2.)/2.0
      !STRK=0.3 ! 0.3 T 
      
      !call MSQUAD twice with 1/2 of the strength and at a distance of +/- alen w.r.t Z0
      !edit: call MSQ twice to call MSQUAD in new coords system
      CALL MSQ(STRK/2.0,APR,X0+alen,x_,y_,z_,Bx,By,Bz)
      CALL MSQ(STRK/2.0,APR,X0-alen,x_,y_,z_,Bx,By,Bz)
      return
      end




