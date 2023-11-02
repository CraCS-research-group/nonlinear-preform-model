!*********************************************************************************************
!************* VUMAT SUBROUTINE FOR A MULTILINEAR ASYMMETRIC CONSTITUTIVE RELATION ***********   
!*************     TO MODEL THE NONLINEAR BENDING BEHAVIOUR OF FABRIC MATERIALS     **********
!                                                                                            *
! @Authors Peter H Broberg, Adam J Thompson, Jonathan P-H Belnoue, Stephen R Hallett,        * 
!           Esben Lindgaard, Brian LV Bak                                                    *                     
! November 2023                                                                              *
! Aalborg University/Universiy of Bristol                                                    *
!                                                                                            *
!**** Notes **********************************************************************************
! This code consisits of the methods described in:                                           *
! P.H. Broberg, E. Lindgaard, A.J. Thompson, J.P-H. Belnoue, S.R. Hallett, and B.L.V. Bak    *
! (2023) "An accurate forming model for capturing the nonlinear material behaviour of        *
! multilayered binder-stabilised fabrics and predicting fibre wrinkling",                    * 
! Manuscript submitted.                                                                      *
!                                                                                            *
! The code is based on the code desribed in the paper:                                       *
! A.J. Thompson, J.P-H. Belnoue, and S.R. Hallett (2020)                                     * 
! "Modelling defect formation in textiles during the double diaphragm forming process",      *
! Composites Part B: Engineering, 202:108357.                                                *
! This code can be found on the repository                                                   *
! https://bristolcompositesinstitute.github.io/HypoDrape                                     *
!                                                                                            *
! The present code contains 2 parts. Part 1 is the implementation of the fibre tracking      *
! algorithm and hypoelastic constitutive law described in Thompson et al. (2020). Part 2 is  *
! the implementation of the asymmetric modulus for describing the nonlinear bending          *
! behaviour described in Broberg et al. (2023). The different parts are highlighted in the   *
! code below.                                                                                *      
!                                                                                            *      
! Please feel free to use and adapt the codes but remember to give proper attribution.       *      
**********************************************************************************************
!**** Input properties ***********************************************************************
! prop(1)-prop(3):   Initial fibre direction 
! prop(4):           Tensile stiffness in direction 1
! prop(5):           Tensile stiffness in direction 2
! prop(6)-props(11): Parameters for the shear stiffness
! prop(12):          Compressive stiffness in direction 1
! prop(13):          Softening stiffness in direction 1
! prop(14):          Softening strain in direction 1
! prop(15):          Locking strain in direction 1
! prop(16):          Compressive stiffness in direction 2
! prop(17):          Softening stiffness in direction 2
! prop(18):          Softening strain in direction 2
! prop(19):          Locking strain in direction 2
! prop(20):          Strain scaling
! prop(21):          Viscous damping
!********************************************************************************************* 
!**** State variables are stored as **********************************************************
! state(*,1)-state(*,4): Stresses in fibre directions
! state(*,5):            Shear angle
! state(*,6):            Strain in fibre direction 1
! state(*,7):            Strain in fibre direction 2
! state(*,8):            Stiffness in direction 1 at previous increment
! state(*,9):            Stiffness in direction 2 at previous increment
!*********************************************************************************************  

      subroutine VUMAT(
! Read only -
     * nblock, ndir, nshr, nstatev, nfieldv, nprops,
     * lanneal, stepTime, totalTime, dt, cmname, coordMp, charLength,
     * props, density, strainInc, relSpinInc,
     * tempOld, stretchOld, defgradOld, fieldOld,
     * stressOld, stateOld, enerInternOld, enerInelasOld,
     * tempNew, stretchNew, defgradNew, fieldNew,
! Write only -
     * stressNew, stateNew, enerInternNew, enerInelasNew)
!
      include 'vaba_param.inc'
!
      character*(*) CMNAME 
      dimension props(nprops), coordMp(nblock,*),
     * charLength(nblock), strainInc(nblock,ndir+nshr),
     * relSpinInc(nblock,nshr), tempOld(nblock),
     * stretchOld(nblock,ndir+nshr),
     * defgradOld(nblock,ndir+nshr+nshr),
     * fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     * stateOld(nblock,nstatev), enerInternOld(nblock),
     * enerInelasOld(nblock), tempNew(nblock),
     * stretchNew(nblock,ndir+nshr),
     * defgradNew(nblock,ndir+nshr+nshr),
     * fieldNew(nblock,nfieldv),
     * stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     * enerInternNew(nblock), enerInelasNew(nblock)	
!
      integer, parameter :: dp = selected_real_kind(p=15, r=307)  
      real(dp), dimension(2,2):: F0, U1, U_1, R1
      real(dp), dimension(2,2):: T1, T2, TT1, TT2, dE, dETT1, dETT2,TdETT1,
     *       TdETT2, sig1, sig2, sig1T, sig2T, sig1TT, sig2TT
      real(dp), dimension(2):: Adir, e01, e02, f01, f02, f1, f2, e1,
     *       e2, f011, f021
      real(dp), dimension(3):: Delta_sig1, Delta_sig2
      real(dp) E11, E22, G12, theta1, theta2, f01_ab, f02_ab, cross1, cross2,
     *       dot1, dot2, Gamma, f1n, f2n, PN1, PN2, PN3, PN4, PN5, PN6, e1n, 
     *       e2n, cos1, cos2, sin1, sin2
      integer, parameter :: ZERO = 0.D0, ONE = 1.0D0
      logical failel	
c      
      do itno=ONE,nblock
!*******************************************************************************
! PART 1 - https://bristolcompositesinstitute.github.io/HypoDrape        
!*******************************************************************************         
!*******************************************************************************
! Grab input parameters of the consititutive law
!*******************************************************************************
! Principle fibre direction 1
      Adir(1) = props(1)
      Adir(2) = props(2)
! modulus along fibre direction 1 and 2
      E11 = props(4)
      E22 = props(5)
! Coeffecients for polynomial defining shear behaviour      
      PN1 = props(6)
      PN2 = props(7)
      PN3 = props(8)
      PN4 = props(9)
      PN5 = props(10)
      PN6 = props(11)
!*******************************************************************************
! Build deformation gradient tensor
      F0 = ZERO
      F0(1,1) = defgradnew(itno,1)
      F0(2,2) = defgradnew(itno,2)
      F0(1,2) = defgradnew(itno,4)	
      F0(2,1) = defgradnew(itno,5)	 

! Build stretch tensor
      U1 = ZERO
      U1(1,1) = stretchNew(itno,1)
      U1(2,2) = stretchNew(itno,2)
      U1(1,2) = stretchNew(itno,4)	
      U1(2,1) = stretchNew(itno,4)	 
      
! Calculate rotation tensor from deformation gradient and stretch tensor
      U_1 = ZERO
      call inverse(U1, U_1, 2)	  
      R1 = matmul(F0, U_1)

! Define initial orthogonal coordinate system of the code work frame
      e01 = (/ONE,ZERO/)
      e02 = (/ZERO,ONE/)
      
! Define initial fibre direction in respect to code workframe
! assumes fibre directions begin orthogonal to one another
      f011(1) = Adir(1)
      f011(2) = Adir(2)
      f021(1) = -Adir(2)
      f021(2) = Adir(1)	  
      
      f1n = sqrt(f011(1)**2+f011(2)**2)
      f2n = sqrt(f021(1)**2+f021(2)**2)

      do i = 1, 2
        f011(i) = f011(i)/f1n
        f021(i) = f021(i)/f2n
      end do

! calculate current orthogonal coordinate system of the code work frame	  
      e1 = matmul(R1, e01)
      e2 = matmul(R1, e02) 

      e1n = sqrt(e1(1)**2+e1(2)**2)
      e2n = sqrt(e2(1)**2+e2(2)**2)

! calculate current fibre directions
      f01 = matmul(F0,f011)
      f02 = matmul(F0,f021)
      
      f01_ab = sqrt(f01(1)**2 + f01(2)**2)
      f02_ab = sqrt(f02(1)**2 + f02(2)**2) 

      f1 = f01/f01_ab
      f2 = f02/f02_ab

! compute angular rotation in fibre direction 1
      dot1 = dot_product(e1, f1)
      cross1  = (e1(1)*f1(2) - e1(2)*f1(1))	
      cos1 = dot1/(e1n*f01_ab)
      sin1 = cross1/(e1n*f01_ab)

      
! compute angular rotation in fibre direction 2      
      dot2 = dot_product(e2, f2)
      cross2  = (e2(1)*f2(2) - e2(2)*f2(1))		
      cos2 = dot2/(e2n*f02_ab)
      sin2 = cross2/(e2n*f02_ab)

! calculate transformation matrix for fibre directions
      T1(1,1) = cos1
      T1(2,2) = cos1
      T1(1,2) = -sin1
      T1(2,1) = sin1	  

      T2(1,1) = cos2
      T2(2,2) = cos2
      T2(1,2) = -sin2
      T2(2,1) = sin2

! store strain increment in tensor
      dE(1,1) = strainInc(itno, 1)
      dE(2,2) = strainInc(itno, 2)
      dE(1,2) = strainInc(itno, 4)
      dE(2,1) = strainInc(itno, 4)  
      
! calculate strain increment in fibre direction 1 and fibre direction 2
      TT1 = transpose(T1)
      TT2 = transpose(T2)

      dETT1 = matmul(TT1, dE)
      dETT2 = matmul(TT2, dE)

      TdETT1 = matmul(dETT1, T1)
      TdETT2 = matmul(dETT2, T2)

! save shear angle, Gamma, as state variable     
      stateNew(itno,5) = stateOld(itno,5) + TdETT1(1,2) + TdETT2(2,1)
      Gamma = abs(stateNew(itno,5))

! Find current shear modulus
      G12 = (PN1*Gamma**5 + PN2*Gamma**4 + PN3*Gamma**3 + 
     *   PN4*Gamma**2 + PN5*Gamma + PN6)
!*******************************************************************************
! PART 2 
!*******************************************************************************       
! accumulate strains in the fibre directions and store as state variables 
      stateNew(itno,6) = stateOld(itno,6) + TdETT1(1,1)
      stateNew(itno,7) = stateOld(itno,7) + TdETT2(2,2)
      
! calcuate the asymmetric modulus
      call asymmetric_modulus(E11, E22, props, stateNew, nProps, 
     1 nBlock, nstatev, itno, stateOld, dt)
      
      call get_regularised_modulus(E11, E22, stateNew, stateOld, 
     * props(21), dt, itno, nblock, nstatev)
      
!*******************************************************************************
! PART 1 - https://bristolcompositesinstitute.github.io/HypoDrape        
!*******************************************************************************           
! compute stress increment in fibre directions
      Delta_sig1(1) = E11 * TdETT1(1,1)
      Delta_sig1(2) = ZERO
      Delta_sig1(3) = G12 * TdETT1(1,2)
    
      Delta_sig2(1) = ZERO
      Delta_sig2(2) = E22 * TdETT2(2,2)
      Delta_sig2(3) = G12 * TdETT2(2,1)

! accumulate stress in the fibre directions and store as state variables
      stateNew(itno,1) = stateOld(itno,1) + Delta_sig1(1)
      stateNew(itno,2) = stateOld(itno,2) + Delta_sig1(3)
      stateNew(itno,3) = stateOld(itno,3) + Delta_sig2(2)
      stateNew(itno,4) = stateOld(itno,4) + Delta_sig2(3) 

! Build stress tensor for each fibre direction with accumulated stress
      sig1(1,1) = stateNew(itno,1)
      sig1(2,2) = ZERO
      sig1(1,2) = stateNew(itno,2)
      sig1(2,1) = stateNew(itno,2)      
      sig2(1,1) = ZERO
      sig2(2,2) = stateNew(itno,3)
      sig2(1,2) = stateNew(itno,4)
      sig2(2,1) = stateNew(itno,4) 

! Transform stress tensors to material frame
      sig1T = matmul(T1, sig1)
      sig2T = matmul(T2, sig2)

      sig1TT = matmul(sig1T,TT1)
      sig2TT = matmul(sig2T,TT2)
        
! Update and return stress values        
      stressNew(itno,1) = sig1TT(1,1) + sig2TT(1,1)
      stressNew(itno,2) = sig1TT(2,2) + sig2TT(2,2)
      stressNew(itno,3) = ZERO
      stressNew(itno,4) = sig1TT(1,2) + sig2TT(1,2)
      
      end do
      return
      end subroutine VUMAT	
!*******************************************************************************
! PART 2
!*******************************************************************************    
      subroutine asymmetric_modulus(E11, E22, props, stateNew, nProps, 
     1 nBlock, nstatev, itno)
      include 'vaba_param.inc'
!============================================================
! Calculate the asymmetric modulus
!-----------------------------------------------------------
! input ...
! props(nProps) - array of properties
! nProps        - dimension
! output ...
! E11           - Computed stiffness in direction 1
! E22           - Computed stiffness in direction 2
! comments ...
! see the description of properties at the code start 
! for details
!===========================================================          
      dimension props(nProps), stateNew(nblock,nstatev)
      
      parameter ( i_prp_E11tens      =  4,
     *            i_prp_E11comp      = 12,
     *            i_prp_E11soft      = 13,
     *            i_prp_eps11soft    = 14,
     *            i_prp_eps11lock    = 15,
     *            i_prp_E22tens      =  5,
     *            i_prp_E22comp      = 16,
     *            i_prp_E22soft      = 17,
     *            i_prp_eps22soft    = 18,
     *            i_prp_eps22lock    = 19,
     *            i_prp_scale_strain = 20)
      
      parameter ( i_sta_eps11 = 6,
     *            i_sta_eps22 = 7)
      
      real*8 E11, E22
! grab the scaling value used for the hyperbolic functions
      scale_strain = props(i_prp_scale_strain)
      
! calculate the modulus in direction 1       
      diffE1_1 = (props(i_prp_E11tens)-props(i_prp_E11comp))/2
      diffE1_2 = (props(i_prp_E11comp)-props(i_prp_E11soft))/2
      diffE1_3 = (props(i_prp_E11soft)-props(i_prp_E11tens))/2
      
      tanh1 = tanh(stateNew(itno,i_sta_eps11) * scale_strain) - 1 
      tanh2 = tanh((stateNew(itno,i_sta_eps11) + props(i_prp_eps11soft))
     *                                          * scale_strain)-1
      tanh3 = tanh((stateNew(itno,i_sta_eps11) + props(i_prp_eps11lock))
     *                                          * scale_strain)-1
      
      E11 = props(i_prp_E11tens) + tanh1 * diffE1_1 + tanh2 * diffE1_2  
     *                                  + tanh3 * diffE1_3 
            
! calculate the modulus in direction 2           
      diffE2_1 = (props(i_prp_E22tens)-props(i_prp_E22comp))/2
      diffE2_2 = (props(i_prp_E22comp)-props(i_prp_E22soft))/2
      diffE2_3 = (props(i_prp_E22soft)-props(i_prp_E22tens))/2
      
      tanh1 = tanh(stateNew(itno,i_sta_eps22) * scale_strain) - 1 
      tanh2 = tanh((stateNew(itno,i_sta_eps22) + props(i_prp_eps22soft))
     *                                          * scale_strain)-1
      tanh3 = tanh((stateNew(itno,i_sta_eps22) + props(i_prp_eps22lock))
     *                                          * scale_strain)-1
      
      E22 = props(i_prp_E22tens) + tanh1 * diffE2_1 + tanh2 * diffE2_2  
     *                                  + tanh3 * diffE2_3 
       
      end subroutine 
      
      subroutine get_regularised_modulus(E11, E22, stateNew, stateOld, 
     * visc_coef, dt, itno, nblock, nstatev)
!============================================================
! Apply damping to the obtained modulus
!-----------------------------------------------------------
! input ...
! E11       - Computed stiffness in direction 1
! E22       - Computed stiffness in direction 2
! visc_coef - Viscous coefficient
! dt        - Size of the time step    
! output ...
! E11       - Regularised stiffness in direction 1
! E22       - Regularised stiffness in direction 2
! comments ...
! the original stiffness E11 and E22 will be destroyed 
! during the calculation
!===========================================================
      include 'vaba_param.inc'
      
      dimension stateNew(nblock,nstatev), stateOld(nblock,nstatev)
      
      real*8 E11, E22
      
      parameter ( i_sta_E11 = 8,
     *            i_sta_E22 = 9)
      
      stateNew(itno,i_sta_E11) = E11 * dt / (visc_coef + dt) +
     *                   stateOld(itno,i_sta_E11) * visc_coef / (visc_coef + dt)
      E11 = stateNew(itno, i_sta_E11)
      
      stateNew(itno,i_sta_E22) = E22 * dt / (visc_coef + dt) +
     *                   stateOld(itno,i_sta_E22) * visc_coef / (visc_coef + dt)
      E22 = stateNew(itno, i_sta_E22)
      
      end subroutine
!*******************************************************************************
! PART 1 - https://bristolcompositesinstitute.github.io/HypoDrape        
!*******************************************************************************          
      subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
      implicit none 
      integer n
      double precision L(n,n), U(n,n), b(n), d(n), x(n), a(n,n), c(n,n)
      double precision coeff
      integer i, j, k 

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
      L=0.0
      U=0.0
      b=0.0

! step 1: forward elimination
      do k=1, n-1
         do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
               a(i,j) = a(i,j)-coeff*a(k,j)
            end do
         end do
      end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
      do i=1,n
        L(i,i) = 1.0
      end do
! U matrix is the upper triangular part of A
      do j=1,n
        do i=1,j
          U(i,j) = a(i,j)
        end do
      end do
c
! Step 3: compute columns of the inverse matrix C
      do k=1,n
        b(k)=1.0
        d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
          d(i)=b(i)
          do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
          end do
        end do
! Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
          x(i) = d(i)
          do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
          end do
          x(i) = x(i)/u(i,i)
        end do
! Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
          c(i,k) = x(i)
        end do
        b(k)=0.0
      end do
      end subroutine inverse
  


