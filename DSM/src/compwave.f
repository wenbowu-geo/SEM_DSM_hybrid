        subroutine comp_excitation_single_force
     &    (maxngrid_r,l,idim_rs_sph,idim_rs_tor,fr,ftheta,fphi,
     &     whole_vector_sph,whole_vector_tor)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing a excitation vector for the given frequency.
c    required subroutines: error_handling
c    required functions: unelastic_factor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
        implicit none
        integer maxngrid_r,ngrid_r,l
        integer idim_rs_sph,idim_rs_tor
        real*8 fr,ftheta,fphi
        complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
        complex*16 whole_vector_tor(maxngrid_r,-2:2)
c other variables
        real*8 b
c constant
        real*8 pi
        parameter ( pi=3.1415926535897932d0)
        if ( l.le.0 ) call error_handling(51)
c **********************************************************************
c Initialing the whole_vector
c **********************************************************************
        call init_complex_array( 10*maxngrid_r,whole_vector_sph )
        call init_complex_array(  5*maxngrid_r,whole_vector_tor )
c **********************************************************************
c Computing the excitation vector
c **********************************************************************
        b=dsqrt(dble(2*l+1)/(4.d0*pi))
        whole_vector_sph(idim_rs_sph,0)
     &    =whole_vector_sph(idim_rs_sph,0)
     &     -b*fr
        whole_vector_sph(idim_rs_sph+1,1)
     &    =whole_vector_sph(idim_rs_sph+1,1)
     &     -0.5d0*b*dcmplx(-ftheta,fphi)
        whole_vector_sph(idim_rs_sph+1,-1)
     &    =whole_vector_sph(idim_rs_sph+1,-1)
     &     -0.5d0*b*dcmplx(ftheta,fphi)
        whole_vector_tor(idim_rs_tor,1)
     &    = whole_vector_tor(idim_rs_tor,1)
     &     -0.5d0*b*dcmplx(fphi,ftheta)
cccccc     &     +0.5d0*b*dcmplx(ftheta,-fphi)
        whole_vector_tor(idim_rs_tor,-1)
     &    = whole_vector_tor(idim_rs_tor,-1)
     &     -0.5d0*b*dcmplx(-fphi,ftheta)
cccc     &     +0.5d0*b*dcmplx(ftheta,-fphi)
ccccccc        print *,"whole_vector",whole_vector_sph(idim_rs_sph,1),
ccccccc     &     whole_vector_sph(idim_rs_sph,-1)
      
        end subroutine comp_excitation_single_force
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine comp_excitation0_single_force
     &    (maxngrid_r,l,idim_rs_sph,fr,whole_vector_sph)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing a excitation vector for the given frequency.
c    required subroutines: error_handling
c    required functions: unelastic_factor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit none
        integer maxngrid_r,l
        integer idim_rs_sph
        real*8 fr
        complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
c other variables
        real*8 b
c constant
        real*8 pi
        parameter ( pi=3.1415926535897932d0)

        if ( l.ne.0 ) call error_handling(52)
c **********************************************************************
c Initialing the whole_vector
c **********************************************************************
        call init_complex_array( 10*maxngrid_r,whole_vector_sph )
c **********************************************************************
c Computing the excitation vector
c **********************************************************************
        b=dsqrt(dble(2*l+1)/(4.d0*pi))
        whole_vector_sph(idim_rs_sph,0)
     &    =whole_vector_sph(idim_rs_sph,0)
     &     -b*fr

        end subroutine comp_excitation0_single_force

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine comp_excitation
     &	  ( maxngrid_r,omega,
     &	    ngrid_r,grid_r,l,source_r,source_mt,
     &	    igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms,
     &	    grid_qkappas,grid_qmus,
     &	    submatrix_I0,submatrix_I1k,submatrix_I1m,
     &	    submatrix_I2,submatrix_I3k,submatrix_I3m,
     &	    submatrix_I4,submatrix_I5k,submatrix_I5m,
     &	    submatrix_I6,submatrix_I7,
     &	    submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &	    idim_rs_sph,idim_rs_tor,
     &	    whole_vector_sph,whole_vector_tor )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing a excitation vector for the given frequency.
c    required subroutines: error_handling
c    required functions: unelastic_factor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
	integer maxngrid_r,ngrid_r,l,igrid_rs
	integer idim_rs_sph,idim_rs_tor
	real*8 grid_r(*)
	real*8 source_r,source_mt(3,3)
	real*8 grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms
	real*8 grid_qkappas,grid_qmus
	real*8 submatrix_I0(4,maxngrid_r)
	real*8 submatrix_I1k(4,maxngrid_r),submatrix_I1m(4,maxngrid_r)
	real*8 submatrix_I2(4,maxngrid_r)
	real*8 submatrix_I3k(4,maxngrid_r),submatrix_I3m(4,maxngrid_r)
	real*8 submatrix_I4(4,maxngrid_r)
	real*8 submatrix_I5k(4,maxngrid_r),submatrix_I5m(4,maxngrid_r)
	real*8 submatrix_I6(4,maxngrid_r)
	real*8 submatrix_I7(4,maxngrid_r)
	real*8 submatrix_I3k_mod(6,maxngrid_r)
	real*8 submatrix_I3m_mod(6,maxngrid_r)
	real*8 submatrix_I4_mod(6,maxngrid_r)
	complex*16 omega
	complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
	complex*16 whole_vector_tor(maxngrid_r,-2:2)
c other variables
	real*8 b1,b2,lsq2,lsq
	complex*16 D1,D2_p,D2_m,D3_p,D3_m,unelastic_factor
	complex*16 factors_qkappa,factors_qmu
c constant
	real*8 pi
	parameter ( pi=3.1415926535897932d0 )
c
	if ( l.le.0 ) call error_handling(51)
c **********************************************************************
c Initialing the whole_vector
c **********************************************************************
	call init_complex_array( 10*maxngrid_r,whole_vector_sph )
	call init_complex_array(  5*maxngrid_r,whole_vector_tor )
c **********************************************************************
c Computing the excitation vector
c **********************************************************************
	b1 = dsqrt( dble(2*l+1) / ( 16.d0 * pi ) )
	b2 = dsqrt( dble(2*l+1)*dble(l-1)*dble(l+2) / ( 64.d0 * pi ) )
	factors_qkappa = unelastic_factor( dble(omega), grid_qkappas )
	factors_qmu = unelastic_factor( dble(omega), grid_qmus )
c        factors_qkappa=cmplx(1.d0,0)
c        factors_qmu = cmplx(1.d0,0)
	lsq2 = dble(l) * dble(l+1)
	lsq = dsqrt( lsq2 )
c --- spheroidal excitations due to the traction diccontinuities
	whole_vector_sph(idim_rs_sph,0)
     &	  = whole_vector_sph(idim_rs_sph,0)
     &	    - b1 * 2.d0
     &	      * ( source_mt(2,2) + source_mt(3,3)
     &	          - 2.d0 * source_mt(1,1)
     &	            * ( grid_Fks * factors_qkappa
     &	                + grid_Fms * factors_qmu )
     &	            / ( grid_Cks * factors_qkappa
     &	                + grid_Cms * factors_qmu )
     &	          ) / source_r
	whole_vector_sph(idim_rs_sph+1,0)
     &	  = whole_vector_sph(idim_rs_sph+1,0)
     &	    - b1 * lsq
     &	      * ( - source_mt(2,2) - source_mt(3,3)
     &	          + 2.d0 * source_mt(1,1)
     &	            * ( grid_Fks * factors_qkappa
     &	                + grid_Fms * factors_qmu )
     &	            / ( grid_Cks * factors_qkappa
     &	                + grid_Cms * factors_qmu )
     &	          ) / source_r
	whole_vector_sph(idim_rs_sph+1,-2)
     &	  = whole_vector_sph(idim_rs_sph+1,-2)
     &	    + b2
     &	      * dcmplx( - source_mt(2,2) + source_mt(3,3),
     &	                - 2.d0 * source_mt(2,3) )
     &	        / source_r
	whole_vector_sph(idim_rs_sph+1,2)
     &	  = whole_vector_sph(idim_rs_sph+1,2)
     &	    + b2
     &	      * dcmplx( - source_mt(2,2) + source_mt(3,3),
     &	                  2.d0 * source_mt(2,3) )
     &	        / source_r
c --- toroidal excitations due to the traction diccontinuities
	whole_vector_tor(idim_rs_tor,-2)
     &	  = whole_vector_tor(idim_rs_tor,-2)
     &	    - b2 * ( dcmplx( - 2.d0 * source_mt(2,3),
     &	                       source_mt(2,2) - source_mt(3,3) ) )
     &	      / source_r
	whole_vector_tor(idim_rs_tor,2)
     &	  = whole_vector_tor(idim_rs_tor,2)
     &	    - b2 * ( dcmplx( - 2.d0 * source_mt(2,3),
     &	                     - source_mt(2,2) + source_mt(3,3) ) )
     &	      / source_r
c --- excitations due to the displacement discontinuities
	D1 = b1 * 2.d0 * source_mt(1,1)
     &	     / ( source_r * source_r 
     &	         * ( grid_Cks * factors_qkappa
     &	             + grid_Cms * factors_qmu )
     &	        )
	D2_p = b1 * dcmplx( -source_mt(1,2), source_mt(1,3) )
     &	     / ( source_r * source_r * grid_Ls * factors_qmu )
	D2_m = b1 * dcmplx( source_mt(1,2), source_mt(1,3) )
     &	     / ( source_r * source_r * grid_Ls * factors_qmu )
	D3_p = b1 * dcmplx( source_mt(1,3), source_mt(1,2) )
     &	     / ( source_r * source_r * grid_Ls * factors_qmu )
	D3_m = b1 * dcmplx( -source_mt(1,3), source_mt(1,2) )
     &	     / ( source_r * source_r * grid_Ls * factors_qmu )
c ---- spheroidal, m=0
	whole_vector_sph(idim_rs_sph,0)
     &	  = whole_vector_sph(idim_rs_sph,0)
     &	      + ( - omega * omega * submatrix_I0(1,igrid_rs)
     &	          + ( submatrix_I1k(1,igrid_rs)
     &	              + 4.d0 * submatrix_I3k(1,igrid_rs)
     &	              + 4.d0 * submatrix_I5k(1,igrid_rs)
     &	            ) * factors_qkappa
     &	          + ( submatrix_I1m(1,igrid_rs)
     &	              + 4.d0 * submatrix_I3m(1,igrid_rs)
     &	              + 4.d0 * submatrix_I5m(1,igrid_rs)
     &	              + lsq2 * submatrix_I6(1,igrid_rs)
     &	              - 4.d0 * submatrix_I7(1,igrid_rs)
     &	            ) * factors_qmu
     &	        ) * D1
	whole_vector_sph(idim_rs_sph+1,0)
     &	  = whole_vector_sph(idim_rs_sph+1,0)
     &	      + ( - lsq * (
     &	              ( submatrix_I3k_mod(1,igrid_rs)
     &	                + 2.d0 * submatrix_I5k(1,igrid_rs)
     &	              ) * factors_qkappa
     &	              + ( submatrix_I3m_mod(1,igrid_rs)
     &	                  - submatrix_I4_mod(1,igrid_rs)
     &	                  + 2.d0 * submatrix_I5m(1,igrid_rs)
     &	                  + submatrix_I6(1,igrid_rs)
     &	                  - 2.d0 * submatrix_I7(1,igrid_rs)
     &	              ) * factors_qmu
     &	        ) ) * D1
	whole_vector_sph(idim_rs_sph+2,0)
     &	  = whole_vector_sph(idim_rs_sph+2,0)
     &	      + ( - omega * omega * submatrix_I0(3,igrid_rs)
     &	          + ( submatrix_I1k(3,igrid_rs)
     &	              + 2.d0 * submatrix_I3k(2,igrid_rs)
     &	              + 2.d0 * submatrix_I3k(3,igrid_rs)
     &	              + 4.d0 * submatrix_I5k(3,igrid_rs)
     &	            ) * factors_qkappa
     &	          + ( submatrix_I1m(3,igrid_rs)
     &	              + 2.d0 * submatrix_I3m(2,igrid_rs)
     &	              + 2.d0 * submatrix_I3m(3,igrid_rs)
     &	              + 4.d0 * submatrix_I5m(3,igrid_rs)
     &	              + lsq2 * submatrix_I6(3,igrid_rs)
     &	              - 4.d0 * submatrix_I7(3,igrid_rs)
     &	            ) * factors_qmu
     &	        ) * D1
	whole_vector_sph(idim_rs_sph+3,0)
     &	  = whole_vector_sph(idim_rs_sph+3,0)
     &	      + ( - lsq * (
     &	              ( submatrix_I3k_mod(3,igrid_rs)
     &	                + 2.d0 * submatrix_I5k(3,igrid_rs)
     &	              ) * factors_qkappa
     &	              + ( submatrix_I3m_mod(3,igrid_rs)
     &	                  - submatrix_I4_mod(3,igrid_rs)
     &	                  + 2.d0 * submatrix_I5m(3,igrid_rs)
     &	                  + submatrix_I6(3,igrid_rs)
     &	                  - 2.d0 * submatrix_I7(3,igrid_rs)
     &	              ) * factors_qmu
     &	        ) ) * D1
c ---- spheroidal, m=-1
	whole_vector_sph(idim_rs_sph,-1)
     &	  = whole_vector_sph(idim_rs_sph,-1)
     &	      + ( - lsq * (
     &	              ( submatrix_I3k_mod(1,igrid_rs)
     &	                + 2.d0 * submatrix_I5k(1,igrid_rs)
     &	              ) * factors_qkappa
     &	              + ( submatrix_I3m_mod(1,igrid_rs)
     &	                  - submatrix_I4_mod(1,igrid_rs)
     &	                  + 2.d0 * submatrix_I5m(1,igrid_rs)
     &	                  + submatrix_I6(1,igrid_rs)
     &	                  - 2.d0 * submatrix_I7(1,igrid_rs)
     &	              ) * factors_qmu
     &	        ) ) * D2_m
	whole_vector_sph(idim_rs_sph+1,-1)
     &	  = whole_vector_sph(idim_rs_sph+1,-1)
     &	      + ( - omega * omega * submatrix_I0(1,igrid_rs)
     &	          + ( lsq2 * submatrix_I5k(1,igrid_rs)
     &	            ) * factors_qkappa
     &	            + ( submatrix_I2(1,igrid_rs)
     &	                - 2.d0 * submatrix_I4(1,igrid_rs)
     &	                + lsq2 * submatrix_I5m(1,igrid_rs)
     &	                + submatrix_I6(1,igrid_rs)
     &	                - 2.d0 * submatrix_I7(1,igrid_rs)
     &	            ) * factors_qmu
     &	        )  * D2_m
	whole_vector_sph(idim_rs_sph+2,-1)
     &	  = whole_vector_sph(idim_rs_sph+2,-1)
     &	      + ( - lsq * (
     &	              ( submatrix_I3k_mod(2,igrid_rs)
     &	                + 2.d0 * submatrix_I5k(3,igrid_rs)
     &	              ) * factors_qkappa
     &	              + ( submatrix_I3m_mod(2,igrid_rs)
     &	                  - submatrix_I4_mod(2,igrid_rs)
     &	                  + 2.d0 * submatrix_I5m(3,igrid_rs)
     &	                  + submatrix_I6(3,igrid_rs)
     &	                  - 2.d0 * submatrix_I7(3,igrid_rs)
     &	              ) * factors_qmu
     &	        ) ) * D2_m
	whole_vector_sph(idim_rs_sph+3,-1)
     &	  = whole_vector_sph(idim_rs_sph+3,-1)
     &	      + ( - omega * omega * submatrix_I0(3,igrid_rs)
     &	          + ( lsq2 * submatrix_I5k(3,igrid_rs)
     &	            ) * factors_qkappa
     &	            + ( submatrix_I2(3,igrid_rs)
     &	                - submatrix_I4(2,igrid_rs)
     &	                - submatrix_I4(3,igrid_rs)
     &	                + lsq2 * submatrix_I5m(3,igrid_rs)
     &	                + submatrix_I6(3,igrid_rs)
     &	                - 2.d0 * submatrix_I7(3,igrid_rs)
     &	            ) * factors_qmu
     &	        )  * D2_m
	whole_vector_sph(idim_rs_sph+4,-1)
     &	  = whole_vector_sph(idim_rs_sph+4,-1)
     &	      + ( - lsq * (
     &	              ( submatrix_I3k_mod(5,igrid_rs)
     &	              ) * factors_qkappa
     &	              + ( submatrix_I3m_mod(5,igrid_rs)
     &	                  - submatrix_I4_mod(5,igrid_rs)
     &	              ) * factors_qmu
     &	        ) ) * D2_m
c ---- spheroidal, m=+1
	whole_vector_sph(idim_rs_sph,1)
     &	  = whole_vector_sph(idim_rs_sph,1)
     &	      + ( - lsq * (
     &	              ( submatrix_I3k_mod(1,igrid_rs)
     &	                + 2.d0 * submatrix_I5k(1,igrid_rs)
     &	              ) * factors_qkappa
     &	              + ( submatrix_I3m_mod(1,igrid_rs)
     &	                  - submatrix_I4_mod(1,igrid_rs)
     &	                  + 2.d0 * submatrix_I5m(1,igrid_rs)
     &	                  + submatrix_I6(1,igrid_rs)
     &	                  - 2.d0 * submatrix_I7(1,igrid_rs)
     &	              ) * factors_qmu
     &	        ) ) * D2_p
	whole_vector_sph(idim_rs_sph+1,1)
     &	  = whole_vector_sph(idim_rs_sph+1,1)
     &	      + ( - omega * omega * submatrix_I0(1,igrid_rs)
     &	          + ( lsq2 * submatrix_I5k(1,igrid_rs)
     &	            ) * factors_qkappa
     &	            + ( submatrix_I2(1,igrid_rs)
     &	                - 2.d0 * submatrix_I4(1,igrid_rs)
     &	                + lsq2 * submatrix_I5m(1,igrid_rs)
     &	                + submatrix_I6(1,igrid_rs)
     &	                - 2.d0 * submatrix_I7(1,igrid_rs)
     &	            ) * factors_qmu
     &	        )  * D2_p
	whole_vector_sph(idim_rs_sph+2,1)
     &	  = whole_vector_sph(idim_rs_sph+2,1)
     &	      + ( - lsq * (
     &	              ( submatrix_I3k_mod(2,igrid_rs)
     &	                + 2.d0 * submatrix_I5k(3,igrid_rs)
     &	              ) * factors_qkappa
     &	              + ( submatrix_I3m_mod(2,igrid_rs)
     &	                  - submatrix_I4_mod(2,igrid_rs)
     &	                  + 2.d0 * submatrix_I5m(3,igrid_rs)
     &	                  + submatrix_I6(3,igrid_rs)
     &	                  - 2.d0 * submatrix_I7(3,igrid_rs)
     &	              ) * factors_qmu
     &	        ) ) * D2_p
	whole_vector_sph(idim_rs_sph+3,1)
     &	  = whole_vector_sph(idim_rs_sph+3,1)
     &	      + ( - omega * omega * submatrix_I0(3,igrid_rs)
     &	          + ( lsq2 * submatrix_I5k(3,igrid_rs)
     &	            ) * factors_qkappa
     &	            + ( submatrix_I2(3,igrid_rs)
     &	                - submatrix_I4(2,igrid_rs)
     &	                - submatrix_I4(3,igrid_rs)
     &	                + lsq2 * submatrix_I5m(3,igrid_rs)
     &	                + submatrix_I6(3,igrid_rs)
     &	                - 2.d0 * submatrix_I7(3,igrid_rs)
     &	            ) * factors_qmu
     &	        )  * D2_p
	whole_vector_sph(idim_rs_sph+4,1)
     &	  = whole_vector_sph(idim_rs_sph+4,1)
     &	      + ( - lsq * (
     &	              ( submatrix_I3k_mod(5,igrid_rs)
     &	              ) * factors_qkappa
     &	              + ( submatrix_I3m_mod(5,igrid_rs)
     &	                  - submatrix_I4_mod(5,igrid_rs)
     &	              ) * factors_qmu
     &	        ) ) * D2_p
c ---- toroidal, m=-1
	whole_vector_tor(idim_rs_tor,-1)
     &	  = whole_vector_tor(idim_rs_tor,-1)
     &	      + ( - omega * omega * submatrix_I0(1,igrid_rs)
     &	          + ( submatrix_I2(1,igrid_rs)
     &	              - 2.d0 * submatrix_I4(1,igrid_rs)
     &	              + submatrix_I6(1,igrid_rs)
     &	              - ( lsq2 - 2.d0 ) * submatrix_I7(1,igrid_rs)
     &	             ) * factors_qmu
     &	         ) * D3_m
	whole_vector_tor(idim_rs_tor+1,-1)
     &	  = whole_vector_tor(idim_rs_tor+1,-1)
     &	      + ( - omega * omega * submatrix_I0(3,igrid_rs)
     &	          + ( submatrix_I2(3,igrid_rs)
     &	              - submatrix_I4(2,igrid_rs)
     &	              - submatrix_I4(3,igrid_rs)
     &	              + submatrix_I6(3,igrid_rs)
     &	              - ( lsq2 - 2.d0 ) * submatrix_I7(3,igrid_rs)
     &	             ) * factors_qmu
     &	         ) * D3_m
c ---- toroidal, m=+1
	whole_vector_tor(idim_rs_tor,1)
     &	  = whole_vector_tor(idim_rs_tor,1)
     &	      + ( - omega * omega * submatrix_I0(1,igrid_rs)
     &	          + ( submatrix_I2(1,igrid_rs)
     &	              - 2.d0 * submatrix_I4(1,igrid_rs)
     &	              + submatrix_I6(1,igrid_rs)
     &	              - ( lsq2 - 2.d0 ) * submatrix_I7(1,igrid_rs)
     &	             ) * factors_qmu
     &	         ) * D3_p
	whole_vector_tor(idim_rs_tor+1,1)
     &	  = whole_vector_tor(idim_rs_tor+1,1)
     &	      + ( - omega * omega * submatrix_I0(3,igrid_rs)
     &	          + ( submatrix_I2(3,igrid_rs)
     &	              - submatrix_I4(2,igrid_rs)
     &	              - submatrix_I4(3,igrid_rs)
     &	              + submatrix_I6(3,igrid_rs)
     &	              - ( lsq2 - 2.d0 ) * submatrix_I7(3,igrid_rs)
     &	             ) * factors_qmu
     &	         ) * D3_p
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine comp_excitation0
     &	  ( maxngrid_r,omega,
     &	    ngrid_r,grid_r,l,source_r,source_mt,
     &	    igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms,
     &	    grid_qkappas,grid_qmus,
     &	    submatrix_I0,submatrix_I1k,submatrix_I1m,
     &	    submatrix_I2,submatrix_I3k,submatrix_I3m,
     &	    submatrix_I4,submatrix_I5k,submatrix_I5m,
     &	    submatrix_I6,submatrix_I7,
     &	    submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &	    idim_rs_sph,idim_rs_tor,
     &	    whole_vector_sph,whole_vector_tor )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing a excitation vector for the given frequency.
c    required subroutines: error_handling
c    required functions: unelastic_factor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
	integer maxngrid_r,ngrid_r,l,igrid_rs
	integer idim_rs_sph,idim_rs_tor
	real*8 grid_r(*)
	real*8 source_r,source_mt(3,3)
	real*8 grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms
	real*8 grid_qkappas,grid_qmus
	real*8 submatrix_I0(4,maxngrid_r)
	real*8 submatrix_I1k(4,maxngrid_r),submatrix_I1m(4,maxngrid_r)
	real*8 submatrix_I2(4,maxngrid_r)
	real*8 submatrix_I3k(4,maxngrid_r),submatrix_I3m(4,maxngrid_r)
	real*8 submatrix_I4(4,maxngrid_r)
	real*8 submatrix_I5k(4,maxngrid_r),submatrix_I5m(4,maxngrid_r)
	real*8 submatrix_I6(4,maxngrid_r)
	real*8 submatrix_I7(4,maxngrid_r)
	real*8 submatrix_I3k_mod(6,maxngrid_r)
	real*8 submatrix_I3m_mod(6,maxngrid_r)
	real*8 submatrix_I4_mod(6,maxngrid_r)
	complex*16 omega
	complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
	complex*16 whole_vector_tor(maxngrid_r,-2:2)
c other variables
	real*8 b1,b2,lsq2
	complex*16 D1,D2_p,D2_m,D3_p,D3_m,unelastic_factor
	complex*16 factors_qkappa,factors_qmu
c constant
	real*8 pi
	parameter ( pi=3.1415926535897932d0 )
c
	if ( l.ne.0 ) call error_handling(52)
c **********************************************************************
c Initialing the whole_vector
c **********************************************************************
	call init_complex_array( 10*maxngrid_r,whole_vector_sph )
	call init_complex_array(  5*maxngrid_r,whole_vector_tor )
c **********************************************************************
c Computing the excitation vector
c **********************************************************************
	b1 = dsqrt( dble(2*l+1) / ( 16.d0 * pi ) )
c	b2 = dsqrt( dble(2*l+1)*dble(l-1)*dble(l+2) / ( 64.d0 * pi ) )
	b2 = 0.d0
c --- excitations due to the traction diccontinuities
	factors_qkappa = unelastic_factor( dble(omega), grid_qkappas )
	factors_qmu = unelastic_factor( dble(omega), grid_qmus )
c        factors_qkappa=cmplx(1.d0,0)
c        factors_qmu = cmplx(1.d0,0)
	lsq2 = dble(l) * dble(l+1)
	whole_vector_sph(idim_rs_sph,0)
     &	  = whole_vector_sph(idim_rs_sph,0)
     &	    - b1 * 2.d0
     &	      * ( source_mt(2,2) + source_mt(3,3)
     &	          - 2.d0 * source_mt(1,1)
     &	            * ( grid_Fks * factors_qkappa
     &	                + grid_Fms * factors_qmu
     &	              )
     &	            / ( grid_Cks * factors_qkappa
     &	                + grid_Cms * factors_qmu
     &	              )
     &	          ) / source_r
c --- excitations due to the displacement discontinuities
	D1 = b1 * 2.d0 * source_mt(1,1)
     &	     / ( source_r * source_r
     &	         * ( grid_Cks * factors_qkappa
     &	             + grid_Cms * factors_qmu ) )
	D2_p = dcmplx( 0.d0 )
	D2_m = dcmplx( 0.d0 )
	D3_p = dcmplx( 0.d0 )
	D3_m = dcmplx( 0.d0 )
c ---- spheroidal, m=0
	whole_vector_sph(idim_rs_sph,0)
     &	  = whole_vector_sph(idim_rs_sph,0)
     &	      + ( - omega * omega * submatrix_I0(1,igrid_rs)
     &	          + ( submatrix_I1k(1,igrid_rs)
     &	              + 4.d0 * submatrix_I3k(1,igrid_rs)
     &	              + 4.d0 * submatrix_I5k(1,igrid_rs)
     &	            ) * factors_qkappa
     &	          + ( submatrix_I1m(1,igrid_rs)
     &	              + 4.d0 * submatrix_I3m(1,igrid_rs)
     &	              + 4.d0 * submatrix_I5m(1,igrid_rs)
     &	              + lsq2 * submatrix_I6(1,igrid_rs)
     &	              - 4.d0 * submatrix_I7(1,igrid_rs)
     &	            ) * factors_qmu
     &	        ) * D1
	whole_vector_sph(idim_rs_sph+1,0)
     &	  = whole_vector_sph(idim_rs_sph+1,0)
     &	      + ( - omega * omega * submatrix_I0(3,igrid_rs)
     &	          + ( submatrix_I1k(3,igrid_rs)
     &	              + 2.d0 * submatrix_I3k(2,igrid_rs)
     &	              + 2.d0 * submatrix_I3k(3,igrid_rs)
     &	              + 4.d0 * submatrix_I5k(3,igrid_rs)
     &	            ) * factors_qkappa
     &	          + ( submatrix_I1m(3,igrid_rs)
     &	              + 2.d0 * submatrix_I3m(2,igrid_rs)
     &	              + 2.d0 * submatrix_I3m(3,igrid_rs)
     &	              + 4.d0 * submatrix_I5m(3,igrid_rs)
     &	              + lsq2 * submatrix_I6(3,igrid_rs)
     &	              - 4.d0 * submatrix_I7(3,igrid_rs)
     &	            ) * factors_qmu
     &	        ) * D1
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine comp_wavefield
     &	  ( maxngrid_r,omega,
     &	    submatrix_I0,submatrix_I1k,submatrix_I1m,
     &	    submatrix_I2,submatrix_I3k,submatrix_I3m,
     &	    submatrix_I4,submatrix_I5k,submatrix_I5m,
     &	    submatrix_I6,submatrix_I7,
     &	    submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &	    ngrid_r,grid_r,grid_mu,grid_qkappa,grid_qmu,l,
     &	    idim1_sph0,idim2_sph,idim1_tor0,idim2_tor,
     &	    idim0,init_npos_sph,init_npos_tor,
     &	    idim_rs_sph,idim_rs_tor,
     &	    idim_station_sph,idim_station_tor,
     &	    whole_matrix_sph,whole_matrix_tor,
     &	    whole_matrix_dr_sph,whole_matrix_dr_tor,
     &	    whole_vector_sph,whole_vector_tor,work_vector,source_type )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing wavefield for the given frequency.
c    required subroutines: init_complex_array
c    required functions: unelastic_factor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
        implicit none
	integer maxngrid_r,ngrid_r,l
	integer idim1_sph0,idim2_sph,idim1_tor0,idim2_tor
	integer idim0,init_npos_sph,init_npos_tor
	integer idim_rs_sph,idim_rs_tor
	integer idim_station_sph,idim_station_tor
	real*8 submatrix_I0(4,maxngrid_r)
	real*8 submatrix_I1k(4,maxngrid_r),submatrix_I1m(4,maxngrid_r)
	real*8 submatrix_I2(4,maxngrid_r)
	real*8 submatrix_I3k(4,maxngrid_r),submatrix_I3m(4,maxngrid_r)
	real*8 submatrix_I4(4,maxngrid_r)
	real*8 submatrix_I5k(4,maxngrid_r),submatrix_I5m(4,maxngrid_r)
	real*8 submatrix_I6(4,maxngrid_r)
	real*8 submatrix_I7(4,maxngrid_r)
	real*8 submatrix_I3k_mod(6,maxngrid_r)
	real*8 submatrix_I3m_mod(6,maxngrid_r)
	real*8 submatrix_I4_mod(6,maxngrid_r)
	real*8 grid_r(*),grid_mu(2,*)
	real*8 grid_qkappa(*),grid_qmu(*)
	complex*16 omega
	complex*16 whole_matrix_sph(4,*),whole_matrix_tor(2,*)
	complex*16 whole_matrix_dr_sph(*),whole_matrix_dr_tor(*)
	complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
	complex*16 whole_vector_tor(maxngrid_r,-2:2)
	complex*16 work_vector(*)
c other variables
	integer idim1_sph,idim1_tor,ir,npos,m,itype_medium
c  WENBO 
        integer source_type,max_m,min_m
c itype_medium=1: solid, itype_medium=0: liquid
	integer ndim_whole_matrix_sph,ndim_whole_matrix_tor,ier
	integer init_grid,end_grid,ns,nq
	real*8 lsq,lsq2,eps
	complex*16 unelastic_factor,factor_qkappa(2),factor_qmu(2)
	data eps / -1.d0 /
c
	if ( l.le.0 ) call error_handling(53)
c **********************************************************************
c Initialing the whole_matrix
c **********************************************************************
	call init_complex_array( 8*maxngrid_r,whole_matrix_sph )
	call init_complex_array( 2*maxngrid_r,whole_matrix_tor )
	call init_complex_array( 2*maxngrid_r,whole_matrix_dr_sph )
	call init_complex_array(   maxngrid_r,whole_matrix_dr_tor )
	idim1_sph = max0( idim0,idim1_sph0 )
	idim1_tor = max0( idim0,idim1_tor0 )
c **********************************************************************
c **********************************************************************
c Spheroidal Component
c **********************************************************************
c **********************************************************************
c **********************************************************************
c constructing the whole_matrix
c **********************************************************************
	lsq2 = dble(l) * dble(l+1)
	lsq  = dsqrt( lsq2 )
c
	factor_qkappa(2)
     &	  = unelastic_factor( dble(omega),grid_qkappa(idim1_sph) )
	factor_qmu(2)
     &	  = unelastic_factor( dble(omega),grid_qmu(idim1_sph) )
	if ( grid_mu(1,idim1_sph)*grid_mu(2,idim1_sph).eq.0.d0 ) then
	  npos = init_npos_sph
	  itype_medium = 0
	  whole_matrix_sph(4,npos)
     &	  = omega * omega / factor_qkappa(2)
     &	      * dcmplx( submatrix_I0(1,idim1_sph)
     &	              )
     &	      - dcmplx( lsq2 * submatrix_I1k(1,idim1_sph)
     &	                  + submatrix_I2(1,idim1_sph)
     &	              )
	else
	  npos = init_npos_sph
	  itype_medium = 1
	  whole_matrix_sph(4,npos)
     &	  = omega * omega
     &	      * dcmplx( submatrix_I0(1,idim1_sph)
     &	              )
     &	      - factor_qkappa(2)
     &	        * dcmplx( submatrix_I1k(1,idim1_sph)
     &	                  + 4.d0 * submatrix_I3k(1,idim1_sph)
     &	                  + 4.d0 * submatrix_I5k(1,idim1_sph)
     &	                )
     &	      - factor_qmu(2)
     &	        * dcmplx( submatrix_I1m(1,idim1_sph)
     &	                  + 4.d0 * submatrix_I3m(1,idim1_sph)
     &	                  + 4.d0 * submatrix_I5m(1,idim1_sph)
     &	                  + lsq2 * submatrix_I6(1,idim1_sph)
     &	                  - 4.d0 * submatrix_I7(1,idim1_sph)
     &	                )
	  whole_matrix_sph(3,npos+1)
     &	  = dcmplx( lsq )
     &	      * (
     &	        factor_qkappa(2)
     &	        * dcmplx( submatrix_I3k_mod(1,idim1_sph)
     &	                  + 2.d0 * submatrix_I5k(1,idim1_sph)
     &	                )
     &	      + factor_qmu(2)
     &	        * dcmplx( submatrix_I3m_mod(1,idim1_sph)
     &	                  - submatrix_I4_mod(1,idim1_sph)
     &	                  + 2.d0 * submatrix_I5m(1,idim1_sph)
     &	                  + submatrix_I6(1,idim1_sph)
     &	                  - 2.d0 * submatrix_I7(1,idim1_sph)
     &	                )
     &	       )
	  whole_matrix_sph(4,npos+1)
     &	  = omega * omega
     &	      * dcmplx( submatrix_I0(1,idim1_sph)
     &	              )
     &	      - factor_qkappa(2)
     &	        * dcmplx( lsq2 * submatrix_I5k(1,idim1_sph)
     &	                )
     &	      - factor_qmu(2)
     &	        * dcmplx( submatrix_I2(1,idim1_sph)
     &	                  - 2.d0 * submatrix_I4(1,idim1_sph)
     &	                  + lsq2 * submatrix_I5m(1,idim1_sph)
     &	                  + submatrix_I6(1,idim1_sph)
     &	                  - 2.d0 * submatrix_I7(1,idim1_sph)
     &	                )
	endif
	do 120 ir=idim1_sph+1,idim2_sph-1
	  factor_qkappa(1)
     &	    = unelastic_factor( dble(omega),grid_qkappa(ir-1) )
	  factor_qmu(1)
     &	    = unelastic_factor( dble(omega),grid_qmu(ir-1) )
	  factor_qkappa(2)
     &	    = unelastic_factor( dble(omega),grid_qkappa(ir) )
	  factor_qmu(2)
     &	    = unelastic_factor( dble(omega),grid_qmu(ir) )
	  if ( grid_mu(1,ir)*grid_mu(2,ir).eq.0.d0 ) then
	    if ( itype_medium.eq.1 ) then
	      npos = npos + 2
	      whole_matrix_sph(2,npos)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(2,ir-1)
     &	                )
     &	        - factor_qkappa(1)
     &	          * dcmplx( submatrix_I1k(2,ir-1)
     &	                + 2.d0 * submatrix_I3k(2,ir-1)
     &	                + 2.d0 * submatrix_I3k(3,ir-1)
     &	                + 4.d0 * submatrix_I5k(2,ir-1)
     &	                  )
     &	        - factor_qmu(1)
     &	          * dcmplx( submatrix_I1m(2,ir-1)
     &	                + 2.d0 * submatrix_I3m(2,ir-1)
     &	                + 2.d0 * submatrix_I3m(3,ir-1)
     &	                + 4.d0 * submatrix_I5m(2,ir-1)
     &	                + lsq2 * submatrix_I6(2,ir-1)
     &	                - 4.d0 * submatrix_I7(2,ir-1)
     &	              )
	      whole_matrix_sph(3,npos)
     &	      = whole_matrix_sph(3,npos)
     &	        + dcmplx( lsq )
     &	        * (
     &	        factor_qkappa(1)
     &	        * dcmplx( submatrix_I3k_mod(2,ir-1)
     &	                  + submatrix_I3k_mod(5,ir-1)
     &	                  + 2.d0 * submatrix_I5k(2,ir-1)
     &	                )
     &	        + factor_qmu(1)
     &	          * dcmplx( submatrix_I3m_mod(2,ir-1)
     &	                    + submatrix_I3m_mod(5,ir-1)
     &	                    - submatrix_I4_mod(2,ir-1)
     &	                    - submatrix_I4_mod(5,ir-1)
     &	                    + 2.d0 * submatrix_I5m(2,ir-1)
     &	                    + submatrix_I6(2,ir-1)
     &	                    - 2.d0 * submatrix_I7(2,ir-1)
     &	                  )
     &	         )
	      whole_matrix_sph(4,npos)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(4,ir-1)
     &	                )
     &	        - factor_qkappa(1)
     &	          * dcmplx( submatrix_I1k(4,ir-1)
     &	                    + 4.d0 * submatrix_I3k(4,ir-1)
     &	                    + 4.d0 * submatrix_I5k(4,ir-1)
     &	                  )
     &	        - factor_qmu(1)
     &	          * dcmplx( submatrix_I1m(4,ir-1)
     &	                    + 4.d0 * submatrix_I3m(4,ir-1)
     &	                    + 4.d0 * submatrix_I5m(4,ir-1)
     &	                    + lsq2 * submatrix_I6(4,ir-1)
     &	                    - 4.d0 * submatrix_I7(4,ir-1)
     &	                   )
	      whole_matrix_sph(1,npos+1)
     &	      = dcmplx( lsq )
     &	        * (
     &	          factor_qkappa(1)
     &	          * dcmplx( submatrix_I3k_mod(3,ir-1)
     &	                    + 2.d0 * submatrix_I5k(2,ir-1)
     &	                  )
     &	          + factor_qmu(1)
     &	            * dcmplx( submatrix_I3m_mod(3,ir-1)
     &	                      - submatrix_I4_mod(3,ir-1)
     &	                      + 2.d0 * submatrix_I5m(2,ir-1)
     &	                      + submatrix_I6(2,ir-1)
     &	                      - 2.d0 * submatrix_I7(2,ir-1)
     &	                    )
     &	           )
	      whole_matrix_sph(2,npos+1)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(2,ir-1)
     &	                )
     &	        - factor_qkappa(1)
     &	          * dcmplx( lsq2 * submatrix_I5k(2,ir-1)
     &	                  )
     &	        - factor_qmu(1)
     &	          * dcmplx( submatrix_I2(2,ir-1)
     &	                    - submatrix_I4(2,ir-1)
     &	                    - submatrix_I4(3,ir-1)
     &	                    + lsq2 * submatrix_I5m(2,ir-1)
     &	                    + submatrix_I6(2,ir-1)
     &	                    - 2.d0 * submatrix_I7(2,ir-1)
     &	                  )
	      whole_matrix_sph(3,npos+1)
     &	      = dcmplx( lsq )
     &	        * (
     &	          factor_qkappa(1)
     &	          * dcmplx( submatrix_I3k_mod(4,ir-1)
     &	                    + submatrix_I3k_mod(6,ir-1)
     &	                    + 2.d0 * submatrix_I5k(4,ir-1)
     &	                  )
     &	          + factor_qmu(1)
     &	            * dcmplx( submatrix_I3m_mod(4,ir-1)
     &	                      + submatrix_I3m_mod(6,ir-1)
     &	                      - submatrix_I4_mod(4,ir-1)
     &	                      - submatrix_I4_mod(6,ir-1)
     &	                      + 2.d0 * submatrix_I5m(4,ir-1)
     &	                      + submatrix_I6(4,ir-1)
     &	                      - 2.d0 * submatrix_I7(4,ir-1)
     &	                    )
     &	          )
	      whole_matrix_sph(4,npos+1)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(4,ir-1)
     &	                )
     &	          - factor_qkappa(1)
     &	            * dcmplx( lsq2 * submatrix_I5k(4,ir-1)
     &	                    )
     &	          - factor_qmu(1)
     &	            * dcmplx( submatrix_I2(4,ir-1)
     &	                      - 2.d0 * submatrix_I4(4,ir-1)
     &	                      + lsq2 * submatrix_I5m(4,ir-1)
     &	                      + submatrix_I6(4,ir-1)
     &	                      - 2.d0 * submatrix_I7(4,ir-1)
     &	                    )
	      npos = npos + 2
	      itype_medium = 0
	      whole_matrix_sph(2,npos)
     &	      = omega * grid_r(ir) * grid_r(ir)
	      whole_matrix_sph(4,npos)
     &	        = omega * omega / factor_qkappa(2)
     &	        * dcmplx( submatrix_I0(1,ir)
     &	                )
     &	        - dcmplx( lsq2 * submatrix_I1k(1,ir)
     &	                  + submatrix_I2(1,ir)
     &	                )
	    else
	      npos = npos + 1
	      itype_medium = 0
	      whole_matrix_sph(3,npos)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(2,ir-1) ) / factor_qkappa(1)
     &	        - dcmplx( lsq2 * submatrix_I1k(2,ir-1)
     &	                  + submatrix_I2(2,ir-1)
     &	                )
	      whole_matrix_sph(4,npos)
     &	        = omega * omega
     &	        * ( dcmplx( submatrix_I0(4,ir-1) ) / factor_qkappa(1)
     &	            + dcmplx( submatrix_I0(1,ir) ) / factor_qkappa(2)
     &	          )
     &	        - dcmplx( lsq2 * submatrix_I1k(4,ir-1)
     &	                  + submatrix_I2(4,ir-1)
     &	                )
     &	        - dcmplx( lsq2 * submatrix_I1k(1,ir)
     &	                  + submatrix_I2(1,ir)
     &	                )
	    endif
	  else
	    if ( itype_medium.eq.0 ) then
	      npos = npos + 1
	      whole_matrix_sph(3,npos)
     &	      = omega * omega / factor_qkappa(1)
     &	        * dcmplx( submatrix_I0(2,ir-1) )
     &	        - dcmplx( lsq2 * submatrix_I1k(2,ir-1)
     &	                  + submatrix_I2(2,ir-1)
     &	                )
	      whole_matrix_sph(4,npos)
     &	      = omega * omega / factor_qkappa(1)
     &	        * dcmplx( submatrix_I0(4,ir-1)
     &	                )
     &	        - dcmplx( lsq2 * submatrix_I1k(4,ir-1)
     &	                  + submatrix_I2(4,ir-1)
     &	                )
	      npos = npos + 1
	      itype_medium = 1
	      whole_matrix_sph(3,npos)
     &	      = - omega * grid_r(ir) * grid_r(ir)
	      whole_matrix_sph(4,npos)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(1,ir)
     &	                )
     &	        - factor_qkappa(2)
     &	        * dcmplx( submatrix_I1k(1,ir)
     &	                  + 4.d0 * submatrix_I3k(1,ir)
     &	                  + 4.d0 * submatrix_I5k(1,ir)
     &	                )
     &	        - factor_qmu(2)
     &	        * dcmplx( submatrix_I1m(1,ir)
     &	                  + 4.d0 * submatrix_I3m(1,ir)
     &	                  + 4.d0 * submatrix_I5m(1,ir)
     &	                  + lsq2 * submatrix_I6(1,ir)
     &	                  - 4.d0 * submatrix_I7(1,ir)
     &	                )
	      whole_matrix_sph(3,npos+1)
     &	      = dcmplx( lsq )
     &	      * (
     &	          factor_qkappa(2)
     &	        * dcmplx( submatrix_I3k_mod(1,ir)
     &	                  + 2.d0 * submatrix_I5k(1,ir)
     &	                )
     &	        + factor_qmu(2)
     &	        * dcmplx( submatrix_I3m_mod(1,ir)
     &	                  - submatrix_I4_mod(1,ir)
     &	                  + 2.d0 * submatrix_I5m(1,ir)
     &	                  + submatrix_I6(1,ir)
     &	                  - 2.d0 * submatrix_I7(1,ir)
     &	                )
     &	       )
	      whole_matrix_sph(4,npos+1)
     &	      = omega * omega
     &	      * dcmplx( submatrix_I0(1,ir)
     &	              )
     &	      - factor_qkappa(2)
     &	        * dcmplx( lsq2 * submatrix_I5k(1,ir)
     &	                )
     &	      - factor_qmu(2)
     &	        * dcmplx( submatrix_I2(1,ir)
     &	                  - 2.d0 * submatrix_I4(1,ir)
     &	                  + lsq2 * submatrix_I5m(1,ir)
     &	                  + submatrix_I6(1,ir)
     &	                  - 2.d0 * submatrix_I7(1,ir)
     &	                )
	    else
	      npos = npos + 2
	      itype_medium = 1
	      whole_matrix_sph(2,npos)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(2,ir-1)
     &	                )
     &	        - factor_qkappa(1)
     &	        * dcmplx( submatrix_I1k(2,ir-1)
     &	                  + 2.d0 * submatrix_I3k(2,ir-1)
     &	                  + 2.d0 * submatrix_I3k(3,ir-1)
     &	                  + 4.d0 * submatrix_I5k(2,ir-1)
     &	                )
     &	        - factor_qmu(1)
     &	        * dcmplx( submatrix_I1m(2,ir-1)
     &	                  + 2.d0 * submatrix_I3m(2,ir-1)
     &	                  + 2.d0 * submatrix_I3m(3,ir-1)
     &	                  + 4.d0 * submatrix_I5m(2,ir-1)
     &	                  + lsq2 * submatrix_I6(2,ir-1)
     &	                  - 4.d0 * submatrix_I7(2,ir-1)
     &	                )
	      whole_matrix_sph(3,npos)
     &	      = whole_matrix_sph(3,npos)
     &	      + dcmplx( lsq )
     &	      * (
     &	        factor_qkappa(1)
     &	        * dcmplx( submatrix_I3k_mod(2,ir-1)
     &	                  + 2.d0 * submatrix_I5k(2,ir-1)
     &	                )
     &	        + factor_qmu(1)
     &	        * dcmplx( submatrix_I3m_mod(2,ir-1)
     &	                  - submatrix_I4_mod(2,ir-1)
     &	                  + 2.d0 * submatrix_I5m(2,ir-1)
     &	                  + submatrix_I6(2,ir-1)
     &	                  - 2.d0 * submatrix_I7(2,ir-1)
     &	                )
     &	       )
	      whole_matrix_sph(4,npos)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(4,ir-1)
     &	                  + submatrix_I0(1,ir)
     &	                )
     &	        - factor_qkappa(1)
     &	        * dcmplx( submatrix_I1k(4,ir-1)
     &	                  + 4.d0 * submatrix_I3k(4,ir-1)
     &	                  + 4.d0 * submatrix_I5k(4,ir-1)
     &	                )
     &	        - factor_qmu(1)
     &	        * dcmplx( submatrix_I1m(4,ir-1)
     &	                  + 4.d0 * submatrix_I3m(4,ir-1)
     &	                  + 4.d0 * submatrix_I5m(4,ir-1)
     &	                  + lsq2 * submatrix_I6(4,ir-1)
     &	                  - 4.d0 * submatrix_I7(4,ir-1)
     &	                )
     &	        - factor_qkappa(2)
     &	        * dcmplx( submatrix_I1k(1,ir)
     &	                  + 4.d0 * submatrix_I3k(1,ir)
     &	                  + 4.d0 * submatrix_I5k(1,ir)
     &	                )
     &	        - factor_qmu(2)
     &	        * dcmplx( submatrix_I1m(1,ir)
     &	                  + 4.d0 * submatrix_I3m(1,ir)
     &	                  + 4.d0 * submatrix_I5m(1,ir)
     &	                  + lsq2 * submatrix_I6(1,ir)
     &	                  - 4.d0 * submatrix_I7(1,ir)
     &	                )
	      whole_matrix_sph(1,npos+1)
     &	      = dcmplx( lsq )
     &	      * (
     &	        factor_qkappa(1)
     &	        * dcmplx( submatrix_I3k_mod(3,ir-1)
     &	                  + 2.d0 * submatrix_I5k(2,ir-1)
     &	                )
     &	        + factor_qmu(1)
     &	        * dcmplx( submatrix_I3m_mod(3,ir-1)
     &	                  - submatrix_I4_mod(3,ir-1)
     &	                  + 2.d0 * submatrix_I5m(2,ir-1)
     &	                  + submatrix_I6(2,ir-1)
     &	                  - 2.d0 * submatrix_I7(2,ir-1)
     &	                )
     &	         )
	      whole_matrix_sph(2,npos+1)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(2,ir-1)
     &	              )
     &	        - factor_qkappa(1)
     &	        * dcmplx( lsq2 * submatrix_I5k(2,ir-1)
     &	                )
     &	        - factor_qmu(1)
     &	        * dcmplx( submatrix_I2(2,ir-1)
     &	                  - submatrix_I4(2,ir-1)
     &	                  - submatrix_I4(3,ir-1)
     &	                  + lsq2 * submatrix_I5m(2,ir-1)
     &	                  + submatrix_I6(2,ir-1)
     &	                  - 2.d0 * submatrix_I7(2,ir-1)
     &	                )
	      whole_matrix_sph(3,npos+1)
     &	      = dcmplx( lsq )
     &	        * (
     &	        factor_qkappa(1)
     &	        * dcmplx( submatrix_I3k_mod(4,ir-1)
     &	                  + 2.d0 * submatrix_I5k(4,ir-1)
     &	                )
     &	        + factor_qmu(1)
     &	        * dcmplx( submatrix_I3m_mod(4,ir-1)
     &	                  - submatrix_I4_mod(4,ir-1)
     &	                  + 2.d0 * submatrix_I5m(4,ir-1)
     &	                  + submatrix_I6(4,ir-1)
     &	                  - 2.d0 * submatrix_I7(4,ir-1)
     &	                )
     &	        + factor_qkappa(2)
     &	        * dcmplx( submatrix_I3k_mod(1,ir)
     &	                  + 2.d0 * submatrix_I5k(1,ir)
     &	                )
     &	        + factor_qmu(2)
     &	        * dcmplx( submatrix_I3m_mod(1,ir)
     &	                  - submatrix_I4_mod(1,ir)
     &	                  + 2.d0 * submatrix_I5m(1,ir)
     &	                  + submatrix_I6(1,ir)
     &	                  - 2.d0 * submatrix_I7(1,ir)
     &	                )
     &	         )
	      whole_matrix_sph(4,npos+1)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(4,ir-1)
     &	                + submatrix_I0(1,ir)
     &	              )
     &	        - factor_qkappa(1)
     &	        * dcmplx( lsq2 * submatrix_I5k(4,ir-1)
     &	                )
     &	        - factor_qmu(1)
     &	        * dcmplx( submatrix_I2(4,ir-1)
     &	                  - 2.d0 * submatrix_I4(4,ir-1)
     &	                  + lsq2 * submatrix_I5m(4,ir-1)
     &	                  + submatrix_I6(4,ir-1)
     &	                  - 2.d0 * submatrix_I7(4,ir-1)
     &	                )
     &	        - factor_qkappa(2)
     &	        * dcmplx( lsq2 * submatrix_I5k(1,ir)
     &	                )
     &	        - factor_qmu(2)
     &	        * dcmplx( submatrix_I2(1,ir)
     &	                  - 2.d0 * submatrix_I4(1,ir)
     &	                  + lsq2 * submatrix_I5m(1,ir)
     &	                  + submatrix_I6(1,ir)
     &	                  - 2.d0 * submatrix_I7(1,ir)
     &	                )
	      whole_matrix_sph(1,npos+2)
     &	      = dcmplx( lsq )
     &	      * (
     &	        factor_qkappa(1)
     &	        * dcmplx( submatrix_I3k_mod(5,ir-1)
     &	                )
     &	        + factor_qmu(1)
     &	        * dcmplx( submatrix_I3m_mod(5,ir-1)
     &	                  - submatrix_I4_mod(5,ir-1)
     &	                )
     &	       )
	      whole_matrix_sph(3,npos+2)
     &	      = dcmplx( lsq )
     &	      * (
     &	        factor_qkappa(1)
     &	        * dcmplx( submatrix_I3k_mod(6,ir-1)
     &	                )
     &	        + factor_qmu(1)
     &	        * dcmplx( submatrix_I3m_mod(6,ir-1)
     &	                  - submatrix_I4_mod(6,ir-1)
     &	                )
     &	       )
	    endif
	  endif
  120	continue
	factor_qkappa(1)
     &	  = unelastic_factor( dble(omega),grid_qkappa(idim2_sph-1) )
	factor_qmu(1)
     &	  = unelastic_factor( dble(omega),grid_qmu(idim2_sph-1) )
	if ( itype_medium.eq.0 ) then
	  npos = npos + 1
	  whole_matrix_sph(3,npos)
     &	    = omega * omega / factor_qkappa(1)
     &	        * dcmplx( submatrix_I0(2,idim2_sph-1)
     &	                )
     &	      - dcmplx( lsq2 * submatrix_I1k(2,idim2_sph-1)
     &	                  + submatrix_I2(2,idim2_sph-1)
     &	              )
	    whole_matrix_sph(4,npos)
     &	      = omega * omega / factor_qkappa(1)
     &	        * dcmplx( submatrix_I0(4,idim2_sph-1)
     &	                )
     &	      - dcmplx( lsq2 * submatrix_I1k(4,idim2_sph-1)
     &	                  + submatrix_I2(4,idim2_sph-1)
     &	                )
	  ndim_whole_matrix_sph = npos
	else
	  npos = npos + 2
	  whole_matrix_sph(2,npos)
     &	  = omega * omega
     &	      * dcmplx( submatrix_I0(2,idim2_sph-1)
     &	              )
     &	    - factor_qkappa(1)
     &	      * dcmplx( submatrix_I1k(2,idim2_sph-1)
     &	                + 2.d0 * submatrix_I3k(2,idim2_sph-1)
     &	                + 2.d0 * submatrix_I3k(3,idim2_sph-1)
     &	                + 4.d0 * submatrix_I5k(2,idim2_sph-1)
     &	              )
     &	    - factor_qmu(1)
     &	      * dcmplx( submatrix_I1m(2,idim2_sph-1)
     &	                + 2.d0 * submatrix_I3m(2,idim2_sph-1)
     &	                + 2.d0 * submatrix_I3m(3,idim2_sph-1)
     &	                + 4.d0 * submatrix_I5m(2,idim2_sph-1)
     &	                + lsq2 * submatrix_I6(2,idim2_sph-1)
     &	                - 4.d0 * submatrix_I7(2,idim2_sph-1)
     &	              )
	  whole_matrix_sph(3,npos)
     &	  = whole_matrix_sph(3,npos)
     &	    + dcmplx( lsq )
     &	    * (
     &	      factor_qkappa(1)
     &	      * dcmplx( submatrix_I3k_mod(2,idim2_sph-1)
     &	                + submatrix_I3k_mod(5,idim2_sph-1)
     &	                + 2.d0 * submatrix_I5k(2,idim2_sph-1)
     &	              )
     &	    + factor_qmu(1)
     &	      * dcmplx( submatrix_I3m_mod(2,idim2_sph-1)
     &	                + submatrix_I3m_mod(5,idim2_sph-1)
     &	                - submatrix_I4_mod(2,idim2_sph-1)
     &	                - submatrix_I4_mod(5,idim2_sph-1)
     &	                + 2.d0 * submatrix_I5m(2,idim2_sph-1)
     &	                + submatrix_I6(2,idim2_sph-1)
     &	                - 2.d0 * submatrix_I7(2,idim2_sph-1)
     &	              )
     &	     )
	  whole_matrix_sph(4,npos)
     &	  = omega * omega
     &	      * dcmplx( submatrix_I0(4,idim2_sph-1)
     &	              )
     &	    - factor_qkappa(1)
     &	      * dcmplx( submatrix_I1k(4,idim2_sph-1)
     &	                + 4.d0 * submatrix_I3k(4,idim2_sph-1)
     &	                + 4.d0 * submatrix_I5k(4,idim2_sph-1)
     &	              )
     &	    - factor_qmu(1)
     &	      * dcmplx( submatrix_I1m(4,idim2_sph-1)
     &	                + 4.d0 * submatrix_I3m(4,idim2_sph-1)
     &	                + 4.d0 * submatrix_I5m(4,idim2_sph-1)
     &	                + lsq2 * submatrix_I6(4,idim2_sph-1)
     &	                - 4.d0 * submatrix_I7(4,idim2_sph-1)
     &	              )
	  whole_matrix_sph(1,npos+1)
     &	  = dcmplx( lsq )
     &	    * (
     &	      factor_qkappa(1)
     &	      * dcmplx( submatrix_I3k_mod(3,idim2_sph-1)
     &	                + 2.d0 * submatrix_I5k(2,idim2_sph-1)
     &	              )
     &	    + factor_qmu(1)
     &	      * dcmplx( submatrix_I3m_mod(3,idim2_sph-1)
     &	                - submatrix_I4_mod(3,idim2_sph-1)
     &	                + 2.d0 * submatrix_I5m(2,idim2_sph-1)
     &	                + submatrix_I6(2,idim2_sph-1)
     &	                - 2.d0 * submatrix_I7(2,idim2_sph-1)
     &	              )
     &	     )
	  whole_matrix_sph(2,npos+1)
     &	  = omega * omega
     &	    * dcmplx( submatrix_I0(2,idim2_sph-1)
     &	            )
     &	    - factor_qkappa(1)
     &	      * dcmplx( lsq2 * submatrix_I5k(2,idim2_sph-1)
     &	              )
     &	    - factor_qmu(1)
     &	      * dcmplx( submatrix_I2(2,idim2_sph-1)
     &	                - submatrix_I4(2,idim2_sph-1)
     &	                - submatrix_I4(3,idim2_sph-1)
     &	                + lsq2 * submatrix_I5m(2,idim2_sph-1)
     &	                + submatrix_I6(2,idim2_sph-1)
     &	                - 2.d0 * submatrix_I7(2,idim2_sph-1)
     &	              )
	  whole_matrix_sph(3,npos+1)
     &	  = dcmplx( lsq )
     &	    * (
     &	      factor_qkappa(1)
     &	      * dcmplx( submatrix_I3k_mod(4,idim2_sph-1)
     &	                + submatrix_I3k_mod(6,idim2_sph-1)
     &	                + 2.d0 * submatrix_I5k(4,idim2_sph-1)
     &	              )
     &	    + factor_qmu(1)
     &	      * dcmplx( submatrix_I3m_mod(4,idim2_sph-1)
     &	                + submatrix_I3m_mod(6,idim2_sph-1)
     &	                - submatrix_I4_mod(4,idim2_sph-1)
     &	                - submatrix_I4_mod(6,idim2_sph-1)
     &	                + 2.d0 * submatrix_I5m(4,idim2_sph-1)
     &	                + submatrix_I6(4,idim2_sph-1)
     &	                - 2.d0 * submatrix_I7(4,idim2_sph-1)
     &	              )
     &	     )
	  whole_matrix_sph(4,npos+1)
     &	  = omega * omega
     &	    * dcmplx( submatrix_I0(4,idim2_sph-1)
     &	            )
     &	    - factor_qkappa(1)
     &	      * dcmplx( lsq2 * submatrix_I5k(4,idim2_sph-1)
     &	              )
     &	    - factor_qmu(1)
     &	      * dcmplx( submatrix_I2(4,idim2_sph-1)
     &	                - 2.d0 * submatrix_I4(4,idim2_sph-1)
     &	                + lsq2 * submatrix_I5m(4,idim2_sph-1)
     &	                + submatrix_I6(4,idim2_sph-1)
     &	                - 2.d0 * submatrix_I7(4,idim2_sph-1)
     &	              )
	  ndim_whole_matrix_sph = npos+1
	endif
c **********************************************************************
c computing the wavefield
c **********************************************************************
c imposing fixed boundary conditions at r=0 and free surface boundary
c conditions at the surface
	if ( ( grid_mu(1,idim1_sph)*grid_mu(2,idim1_sph).ne.0.d0 )
     &	     .and.( grid_r(idim1_sph).eq.0.d0 ) ) then
	  init_grid = max0(init_npos_sph,3)
	else
	  init_grid = init_npos_sph
	endif
	if ( grid_mu(1,idim2_sph-1)*grid_mu(2,idim2_sph-1).eq.0.d0 )
     &	  then
	  end_grid = ndim_whole_matrix_sph - 1
	else
	  end_grid = ndim_whole_matrix_sph
	endif
c  WENBO
        if(source_type.eq.1) then
          max_m=2
          min_m=-2
        else if(source_type.eq.2) then
          max_m=1
          min_m=-1
        end if
c  WENBO
	m = max0(-l,min_m)
	ns = idim_rs_sph - init_grid + 1
c	if ( mod(l,100).eq.0 ) then
c	  nq = end_grid - init_grid + 1
c	else
cWENBO
c
c	  nq = min0( end_grid - idim_station_sph + 1,
c     &	             end_grid - init_grid + 1 )

!         nq = min0( end_grid - idim_station_sph + 1 + 2,
!     &              end_grid - init_grid + 1 )

c	endif
!solving all the station grids, even the ones whose amplitudes are below 
!    the threshold (parameter eps in subroutine check_amp_significance)
!         nq=max0(end_grid-idim_station_sph+1, end_grid-init_grid+1+2)
         nq=end_grid-init_grid+1+2

	call dclisb(whole_matrix_sph(1,init_grid),
     &	            end_grid - init_grid + 1,
     &	            3,4,ns,nq,whole_vector_sph(init_grid,m),
     &	            eps,whole_matrix_dr_sph(init_grid),
     &	            work_vector(init_grid),ier)
c WENBO
	do 200 m=max0(-l,min_m)+1,min0(l,max_m)
	  call dcsbsub(whole_matrix_sph(1,init_grid),
     &	               end_grid - init_grid + 1,
     &	               3,4,ns,nq,whole_vector_sph(init_grid,m),
     &	               eps,whole_matrix_dr_sph(init_grid),
     &	               work_vector(init_grid),ier)
  200	continue
c
c **********************************************************************
c **********************************************************************
c Toroidal Component
c **********************************************************************
c **********************************************************************
c **********************************************************************
c constructing the whole_matrix
c **********************************************************************
	factor_qmu(2)
     &	  = unelastic_factor( dble(omega),grid_qmu(idim1_tor) )
	npos = init_npos_tor
	whole_matrix_tor(2,npos)
     &	  = omega * omega
     &	      * dcmplx( submatrix_I0(1,idim1_tor)
     &	              )
     &	      - factor_qmu(2)
     &	        * dcmplx( submatrix_I2(1,idim1_tor)
     &	                  - 2.d0 * submatrix_I4(1,idim1_tor)
     &	                  + submatrix_I6(1,idim1_tor)
     &	                  + ( lsq2 - 2.d0 ) * submatrix_I7(1,idim1_tor)
     &	                )
	do 320 ir=idim1_tor+1,idim2_tor-1
	  factor_qmu(1)
     &	    = unelastic_factor( dble(omega),grid_qmu(ir-1) )
	  factor_qmu(2)
     &	    = unelastic_factor( dble(omega),grid_qmu(ir) )
	  npos = npos + 1
	  whole_matrix_tor(1,npos)
     &	    = omega * omega
     &	        * dcmplx(  submatrix_I0(2,ir-1)
     &	                )
     &	        - factor_qmu(1)
     &	          * dcmplx( submatrix_I2(2,ir-1)
     &	                      - submatrix_I4(2,ir-1)
     &	                      - submatrix_I4(3,ir-1)
     &	                      + submatrix_I6(2,ir-1)
     &	                      + ( lsq2 - 2.d0 )
     &	                        * submatrix_I7(2,ir-1)
     &	                  )
	  whole_matrix_tor(2,npos)
     &	    = omega * omega
     &	        * dcmplx(   submatrix_I0(4,ir-1)
     &	                  + submatrix_I0(1,ir)
     &	                )
     &	        - factor_qmu(1)
     &	          * dcmplx( submatrix_I2(4,ir-1)
     &	                      - 2.d0 * submatrix_I4(4,ir-1)
     &	                      + submatrix_I6(4,ir-1)
     &	                      + ( lsq2 - 2.d0 )
     &	                        * submatrix_I7(4,ir-1)
     &	                  )
     &	        - factor_qmu(2)
     &	          * dcmplx( submatrix_I2(1,ir)
     &	                      - 2.d0 * submatrix_I4(1,ir)
     &	                      + submatrix_I6(1,ir)
     &	                      + ( lsq2 - 2.d0 )
     &	                        * submatrix_I7(1,ir)
     &	                  )
  320	continue
	factor_qmu(1)
     &	  = unelastic_factor( dble(omega),
     &	                      grid_qmu(idim2_tor-1) )
	npos = npos + 1
	whole_matrix_tor(1,npos)
     &	  = omega * omega
     &	      * dcmplx(   submatrix_I0(2,idim2_tor-1)
     &	              )
     &	      - factor_qmu(1)
     &	        * dcmplx( submatrix_I2(2,idim2_tor-1)
     &	                  - submatrix_I4(2,idim2_tor-1)
     &	                  - submatrix_I4(3,idim2_tor-1)
     &	                  + submatrix_I6(2,idim2_tor-1)
     &	                  + ( lsq2 - 2.d0 )
     &	                    * submatrix_I7(2,idim2_tor-1) )
	whole_matrix_tor(2,npos)
     &	  = omega * omega
     &	      * dcmplx(  submatrix_I0(4,idim2_tor-1)
     &	              )
     &	      - factor_qmu(1)
     &	        * dcmplx( submatrix_I2(4,idim2_tor-1)
     &	                  - 2.d0 * submatrix_I4(4,idim2_tor-1)
     &	                  + submatrix_I6(4,idim2_tor-1)
     &	                  + ( lsq2 - 2.d0 )
     &	                    * submatrix_I7(4,idim2_tor-1) )
	ndim_whole_matrix_tor = npos
c **********************************************************************
c computing the wavefield
c **********************************************************************
c WENBO
	m = max0(-l,min_m)
	init_grid = init_npos_tor
	end_grid = ndim_whole_matrix_tor
	ns = idim_rs_tor - init_grid + 1
c	if ( mod(l,100).eq.0 ) then
c	  nq = end_grid - init_grid + 1
c	else
cWENBO
c	  nq = min0( end_grid - idim_station_tor + 1,
c     &	             end_grid - init_grid + 1 )
!         nq = min0( end_grid - idim_station_tor + 1 +2,
!     &              end_grid - init_grid + 1 )

c	endif
!solving all the station grids, even the ones whose amplitudes are below 
!    the threshold (parameter eps in subroutine check_amp_significance)
!         nq=max0(end_grid - init_grid + 1,end_grid - init_grid + 1+2)
         nq=end_grid  - init_grid + 1+2

	call dclisb(whole_matrix_tor(1,init_grid),
     &	            end_grid - init_grid + 1,
     &	            1,2,ns,nq,whole_vector_tor(init_grid,m),
     &	            eps,whole_matrix_dr_tor(init_grid),
     &	            work_vector(init_grid),ier)
c WENBO
	do 400 m=max0(-l,min_m)+1,min0(l,max_m)
	  call dcsbsub(whole_matrix_tor(1,init_grid),
     &	               end_grid - init_grid + 1,
     &	               1,2,ns,nq,whole_vector_tor(init_grid,m),
     &	               eps,whole_matrix_dr_tor(init_grid),
     &	               work_vector(init_grid),ier)
  400	continue
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine comp_wavefield0
     &	  ( maxngrid_r,omega,
     &	    submatrix_I0,submatrix_I1k,submatrix_I1m,
     &	    submatrix_I2,submatrix_I3k,submatrix_I3m,
     &	    submatrix_I4,submatrix_I5k,submatrix_I5m,
     &	    submatrix_I6,submatrix_I7,
     &	    submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &	    ngrid_r,grid_r,grid_mu,grid_qkappa,grid_qmu,l,
     &	    idim1_sph,idim2_sph,idim1_tor,idim2_tor,
     &	    idim0,init_npos_sph,init_npos_tor,
     &	    idim_rs_sph,idim_rs_tor,
     &	    idim_station_sph,idim_station_tor,
     &	    whole_matrix_sph,whole_matrix_tor,
     &	    whole_matrix_dr_sph,whole_matrix_dr_tor,
     &	    whole_vector_sph,whole_vector_tor,work_vector )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing wavefield for the given frequency.
c    required subroutines: init_complex_array,error_handling
c    required functions: unelastic_factor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
        implicit none
	integer maxngrid_r,ngrid_r,l
	integer idim1_sph,idim2_sph,idim1_tor,idim2_tor
	integer idim0,init_npos_sph,init_npos_tor
	integer idim_rs_sph,idim_rs_tor
	integer idim_station_sph,idim_station_tor
	real*8 submatrix_I0(4,maxngrid_r)
	real*8 submatrix_I1k(4,maxngrid_r),submatrix_I1m(4,maxngrid_r)
	real*8 submatrix_I2(4,maxngrid_r)
	real*8 submatrix_I3k(4,maxngrid_r),submatrix_I3m(4,maxngrid_r)
	real*8 submatrix_I4(4,maxngrid_r)
	real*8 submatrix_I5k(4,maxngrid_r),submatrix_I5m(4,maxngrid_r)
	real*8 submatrix_I6(4,maxngrid_r)
	real*8 submatrix_I7(4,maxngrid_r)
	real*8 submatrix_I3k_mod(6,maxngrid_r)
	real*8 submatrix_I3m_mod(6,maxngrid_r)
	real*8 submatrix_I4_mod(6,maxngrid_r)
	real*8 grid_r(*),grid_mu(2,*)
	real*8 grid_qkappa(*),grid_qmu(*)
	complex*16 omega
	complex*16 whole_matrix_sph(4,*),whole_matrix_tor(2,*)
	complex*16 whole_matrix_dr_sph(*),whole_matrix_dr_tor(*)
	complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
	complex*16 whole_vector_tor(maxngrid_r,-2:2)
	complex*16 work_vector(*)
c other variables
	integer ir,npos,m,itype_medium
c itype_medium=1: solid, itype_medium=0: liquid
	integer ndim_whole_matrix_sph,ndim_whole_matrix_tor,ier
	integer init_grid,end_grid,ns,nq
	real*8 lsq,lsq2,eps
	complex*16 unelastic_factor,factor_qkappa(2),factor_qmu(2)
	data eps / -1.d0 /
c
	if ( l.ne.0 ) call error_handling(54)
c **********************************************************************
c Initialing the whole_matrix
c **********************************************************************
	call init_complex_array( 8*maxngrid_r,whole_matrix_sph )
	call init_complex_array( 2*maxngrid_r,whole_matrix_dr_sph )
	idim0 = 1
	init_npos_sph = 1
	init_npos_tor = 1
c **********************************************************************
c **********************************************************************
c Spheroidal Component
c **********************************************************************
c **********************************************************************
c **********************************************************************
c constructing the whole_matrix
c **********************************************************************
	lsq2 = dble(l) * dble(l+1)
	lsq  = dsqrt( lsq2 )
c
	factor_qkappa(2)
     &	  = unelastic_factor( dble(omega),grid_qkappa(idim1_sph) )
	factor_qmu(2)
     &	  = unelastic_factor( dble(omega),grid_qmu(idim1_sph) )
	if ( grid_mu(1,idim1_sph)*grid_mu(2,idim1_sph).eq.0.d0 ) then
	  npos = 1
	  itype_medium = 0
	  whole_matrix_sph(2,npos)
     &	  = omega * omega / factor_qkappa(2)
     &	      * dcmplx( submatrix_I0(1,idim1_sph)
     &	              )
     &	    - dcmplx( lsq2 * submatrix_I1k(1,idim1_sph)
     &	                  + submatrix_I2(1,idim1_sph)
     &	            )
	else
	  npos = 1
	  itype_medium = 1
	  whole_matrix_sph(2,npos)
     &	  = omega * omega
     &	      * dcmplx( submatrix_I0(1,idim1_sph)
     &	              )
     &	      - factor_qkappa(2)
     &	        * dcmplx( submatrix_I1k(1,idim1_sph)
     &	                  + 4.d0 * submatrix_I3k(1,idim1_sph)
     &	                  + 4.d0 * submatrix_I5k(1,idim1_sph)
     &	                )
     &	      - factor_qmu(2)
     &	        * dcmplx( submatrix_I1m(1,idim1_sph)
     &	                  + 4.d0 * submatrix_I3m(1,idim1_sph)
     &	                  + 4.d0 * submatrix_I5m(1,idim1_sph)
     &	                  + lsq2 * submatrix_I6(1,idim1_sph)
     &	                  - 4.d0 * submatrix_I7(1,idim1_sph)
     &	                )
	endif
	do 120 ir=idim1_sph+1,idim2_sph-1
	  factor_qkappa(1)
     &	      = unelastic_factor( dble(omega),grid_qkappa(ir-1) )
	  factor_qmu(1)
     &	      = unelastic_factor( dble(omega),grid_qmu(ir-1) )
	  factor_qkappa(2)
     &	    = unelastic_factor( dble(omega),grid_qkappa(ir) )
	  factor_qmu(2)
     &	    = unelastic_factor( dble(omega),grid_qmu(ir) )
	  if ( grid_mu(1,ir)*grid_mu(2,ir).eq.0.d0 ) then
	    if ( itype_medium.eq.1 ) then
	      npos = npos + 1
	      whole_matrix_sph(1,npos)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(2,ir-1)
     &	                )
     &	        - factor_qkappa(1)
     &	          * dcmplx( submatrix_I1k(2,ir-1)
     &	                + 2.d0 * submatrix_I3k(2,ir-1)
     &	                + 2.d0 * submatrix_I3k(3,ir-1)
     &	                + 4.d0 * submatrix_I5k(2,ir-1)
     &	                  )
     &	        - factor_qmu(1)
     &	          * dcmplx( submatrix_I1m(2,ir-1)
     &	                + 2.d0 * submatrix_I3m(2,ir-1)
     &	                + 2.d0 * submatrix_I3m(3,ir-1)
     &	                + 4.d0 * submatrix_I5m(2,ir-1)
     &	                + lsq2 * submatrix_I6(2,ir-1)
     &	                - 4.d0 * submatrix_I7(2,ir-1)
     &	              )
	      whole_matrix_sph(2,npos)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(4,ir-1)
     &	                )
     &	        - factor_qkappa(1)
     &	          * dcmplx( submatrix_I1k(4,ir-1)
     &	                    + 4.d0 * submatrix_I3k(4,ir-1)
     &	                    + 4.d0 * submatrix_I5k(4,ir-1)
     &	                  )
     &	        - factor_qmu(1)
     &	          * dcmplx( submatrix_I1m(4,ir-1)
     &	                    + 4.d0 * submatrix_I3m(4,ir-1)
     &	                    + 4.d0 * submatrix_I5m(4,ir-1)
     &	                    + lsq2 * submatrix_I6(4,ir-1)
     &	                    - 4.d0 * submatrix_I7(4,ir-1)
     &	                   )
	      npos = npos + 1
	      itype_medium = 0
	      whole_matrix_sph(1,npos)
     &	      = omega * grid_r(ir) * grid_r(ir)
	      whole_matrix_sph(2,npos)
     &	      = omega * omega
     &	        * (
     &	          dcmplx( submatrix_I0(1,ir)   ) / factor_qkappa(2)
     &	          )
     &	      - dcmplx( lsq2 * submatrix_I1k(1,ir)
     &	                  + submatrix_I2(1,ir)
     &	              )
	    else
	      npos = npos + 1
	      whole_matrix_sph(1,npos)
     &	      = omega * omega / factor_qkappa(1)
     &	        * dcmplx( submatrix_I0(2,ir-1) )
     &	        - dcmplx( lsq2 * submatrix_I1k(2,ir-1)
     &	                  + submatrix_I2(2,ir-1)
     &	                )
	      whole_matrix_sph(2,npos)
     &	      = omega * omega
     &	        * (
     &	          dcmplx( submatrix_I0(4,ir-1) ) / factor_qkappa(1)
     &	          + dcmplx( submatrix_I0(1,ir)   ) / factor_qkappa(2)
     &	          )
     &	      - dcmplx( lsq2 * submatrix_I1k(4,ir-1)
     &	                  + submatrix_I2(4,ir-1)
     &	              )
     &	      - dcmplx( lsq2 * submatrix_I1k(1,ir)
     &	                  + submatrix_I2(1,ir)
     &	              )
	    endif
	  else
	    if ( itype_medium.eq.0 ) then
	      npos = npos + 1
	      whole_matrix_sph(1,npos)
     &	      = omega * omega / factor_qkappa(1)
     &	        * dcmplx( submatrix_I0(2,ir-1)
     &	                )
     &	        - dcmplx( lsq2 * submatrix_I1k(2,ir-1)
     &	                  + submatrix_I2(2,ir-1)
     &	                )
	      whole_matrix_sph(2,npos)
     &	      = omega * omega / factor_qkappa(1)
     &	        * dcmplx( submatrix_I0(4,ir-1)
     &	                )
     &	      - dcmplx( lsq2 * submatrix_I1k(4,ir-1)
     &	                  + submatrix_I2(4,ir-1)
     &	              )
	      npos = npos + 1
	      itype_medium = 1
	      whole_matrix_sph(1,npos)
     &	      = - omega * grid_r(ir) * grid_r(ir)
	      whole_matrix_sph(2,npos)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(1,ir)
     &	                )
     &	        - factor_qkappa(2)
     &	        * dcmplx( submatrix_I1k(1,ir)
     &	                  + 4.d0 * submatrix_I3k(1,ir)
     &	                  + 4.d0 * submatrix_I5k(1,ir)
     &	                )
     &	        - factor_qmu(2)
     &	        * dcmplx( submatrix_I1m(1,ir)
     &	                  + 4.d0 * submatrix_I3m(1,ir)
     &	                  + 4.d0 * submatrix_I5m(1,ir)
     &	                  + lsq2 * submatrix_I6(1,ir)
     &	                  - 4.d0 * submatrix_I7(1,ir)
     &	                )
	    else
	      npos = npos + 1
	      itype_medium = 1
	      whole_matrix_sph(1,npos)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(2,ir-1)
     &	                )
     &	        - factor_qkappa(1)
     &	        * dcmplx( submatrix_I1k(2,ir-1)
     &	                  + 2.d0 * submatrix_I3k(2,ir-1)
     &	                  + 2.d0 * submatrix_I3k(3,ir-1)
     &	                  + 4.d0 * submatrix_I5k(2,ir-1)
     &	                )
     &	        - factor_qmu(1)
     &	        * dcmplx( submatrix_I1m(2,ir-1)
     &	                  + 2.d0 * submatrix_I3m(2,ir-1)
     &	                  + 2.d0 * submatrix_I3m(3,ir-1)
     &	                  + 4.d0 * submatrix_I5m(2,ir-1)
     &	                  + lsq2 * submatrix_I6(2,ir-1)
     &	                  - 4.d0 * submatrix_I7(2,ir-1)
     &	                )
	      whole_matrix_sph(2,npos)
     &	      = omega * omega
     &	        * dcmplx( submatrix_I0(4,ir-1)
     &	                  + submatrix_I0(1,ir)
     &	                )
     &	        - factor_qkappa(1)
     &	        * dcmplx( submatrix_I1k(4,ir-1)
     &	                  + 4.d0 * submatrix_I3k(4,ir-1)
     &	                  + 4.d0 * submatrix_I5k(4,ir-1)
     &	                )
     &	        - factor_qmu(1)
     &	        * dcmplx( submatrix_I1m(4,ir-1)
     &	                  + 4.d0 * submatrix_I3m(4,ir-1)
     &	                  + 4.d0 * submatrix_I5m(4,ir-1)
     &	                  + lsq2 * submatrix_I6(4,ir-1)
     &	                  - 4.d0 * submatrix_I7(4,ir-1)
     &	                )
     &	        - factor_qkappa(2)
     &	        * dcmplx( submatrix_I1k(1,ir)
     &	                  + 4.d0 * submatrix_I3k(1,ir)
     &	                  + 4.d0 * submatrix_I5k(1,ir)
     &	                )
     &	        - factor_qmu(2)
     &	        * dcmplx( submatrix_I1m(1,ir)
     &	                  + 4.d0 * submatrix_I3m(1,ir)
     &	                  + 4.d0 * submatrix_I5m(1,ir)
     &	                  + lsq2 * submatrix_I6(1,ir)
     &	                  - 4.d0 * submatrix_I7(1,ir)
     &	                )
	    endif
	  endif
  120	continue
	factor_qkappa(1)
     &	  = unelastic_factor( dble(omega),grid_qkappa(idim2_sph-1) )
	factor_qmu(1)
     &	  = unelastic_factor( dble(omega),grid_qmu(idim2_sph-1) )
	if ( itype_medium.eq.0 ) then
	  npos = npos + 1
	  whole_matrix_sph(1,npos)
     &	    = omega * omega / factor_qkappa(1)
     &	        * dcmplx( submatrix_I0(2,idim2_sph-1)
     &	                )
     &	      - dcmplx( lsq2 * submatrix_I1k(2,idim2_sph-1)
     &	                  + submatrix_I2(2,idim2_sph-1)
     &	              )
	  whole_matrix_sph(2,npos)
     &	      = omega * omega / factor_qkappa(1)
     &	        * dcmplx( submatrix_I0(4,idim2_sph-1)
     &	                )
     &	      - dcmplx( lsq2 * submatrix_I1k(4,idim2_sph-1)
     &	                + submatrix_I2(4,idim2_sph-1)
     &	              )
	  ndim_whole_matrix_sph = npos
	else
	  npos = npos + 1
	  whole_matrix_sph(1,npos)
     &	  = omega * omega
     &	      * dcmplx( submatrix_I0(2,idim2_sph-1)
     &	              )
     &	    - factor_qkappa(1)
     &	      * dcmplx( submatrix_I1k(2,idim2_sph-1)
     &	                + 2.d0 * submatrix_I3k(2,idim2_sph-1)
     &	                + 2.d0 * submatrix_I3k(3,idim2_sph-1)
     &	                + 4.d0 * submatrix_I5k(2,idim2_sph-1)
     &	              )
     &	    - factor_qmu(1)
     &	      * dcmplx( submatrix_I1m(2,idim2_sph-1)
     &	                + 2.d0 * submatrix_I3m(2,idim2_sph-1)
     &	                + 2.d0 * submatrix_I3m(3,idim2_sph-1)
     &	                + 4.d0 * submatrix_I5m(2,idim2_sph-1)
     &	                + lsq2 * submatrix_I6(2,idim2_sph-1)
     &	                - 4.d0 * submatrix_I7(2,idim2_sph-1)
     &	              )
	  whole_matrix_sph(2,npos)
     &	  = omega * omega
     &	      * dcmplx( submatrix_I0(4,idim2_sph-1)
     &	              )
     &	    - factor_qkappa(1)
     &	      * dcmplx( submatrix_I1k(4,idim2_sph-1)
     &	                + 4.d0 * submatrix_I3k(4,idim2_sph-1)
     &	                + 4.d0 * submatrix_I5k(4,idim2_sph-1)
     &	              )
     &	    - factor_qmu(1)
     &	      * dcmplx( submatrix_I1m(4,idim2_sph-1)
     &	                + 4.d0 * submatrix_I3m(4,idim2_sph-1)
     &	                + 4.d0 * submatrix_I5m(4,idim2_sph-1)
     &	                + lsq2 * submatrix_I6(4,idim2_sph-1)
     &	                - 4.d0 * submatrix_I7(4,idim2_sph-1)
     &	              )
	  ndim_whole_matrix_sph = npos
	endif
c **********************************************************************
c computing the wavefield
c **********************************************************************
c imposing fixed boundary conditions at r=0 and free surface boundary
c conditions at the surface
	if ( ( grid_mu(1,idim1_sph)*grid_mu(2,idim1_sph).ne.0.d0 )
     &	     .and.( grid_r(idim1_sph).eq.0.d0 ) ) then
	  init_grid = 2
	else
	  init_grid = 1
	endif
	if ( grid_mu(1,idim2_sph-1)*grid_mu(2,idim2_sph-1).eq.0.d0 )
     &	  then
	  end_grid = ndim_whole_matrix_sph - 1
	else
	  end_grid = ndim_whole_matrix_sph
	endif
	m = 0
	ns = idim_rs_sph - init_grid + 1
	if ( mod(l,100).eq.0 ) then
	  nq = end_grid - init_grid + 1
	else
	  nq = end_grid - idim_station_sph + 1
	endif
c WENBO

	call dclisb(whole_matrix_sph(1,init_grid),
     &	            end_grid - init_grid + 1,
     &	            1,4,ns,nq,whole_vector_sph(init_grid,m),
     &	            eps,whole_matrix_dr_sph(init_grid),
     &	            work_vector(init_grid),ier)

c
	ndim_whole_matrix_tor = 0
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	complex*16 function unelastic_factor( omega,qmu )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing the unelastic factor (ratio between complex mu and
c real mu) for the given quality factor, qmu.
c    required subroutines: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c input/output variables
	real*8 omega,qmu
c other variables
	real*8 vr,vi
c constants
	real*8 pi
	parameter ( pi=3.1415926535897932d0 )
c
c **********************************************************************
c computing the complex velocity for the given qmu
c **********************************************************************
	if ( (omega.eq.0.d0).or.(qmu.lt.0.d0) ) then
	  vr = 1.d0
	  vi = 0.d0
	else
	  vr = 1.d0
     &	       + dlog( omega / ( 2.d0 * pi ) ) / ( pi * qmu )
	  vi = 1.d0 / ( 2.d0 * qmu )
	endif
c
c **********************************************************************
c computing the unelastic factor
c **********************************************************************
	unelastic_factor = dcmplx( vr, vi ) * dcmplx( vr, vi )
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine dclisb(a, n, nud, n1, np, nq, b, eps, dr, z, ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
*  simultaneous linear equations with real symmetric positive definite *
*      band matrix by cholesky method.                                 *
*  parameters                                                          *
*    (1) a : 2-dim. array containing the matrix.                       *
*    (2) n : order of the matrix.                                      *
*    (3) nud : size of band's half width.                              *
*    (4) n1 : row size of the array a in the 'dimension' statement.    *
*    (5) b : 1-dim. array containing the right hand side vector.       *
*    (6) eps : parameter to check singurarity off the matrix           *
*              standard value = 1.0d-14                                *
*    (7) dr : 1-dim. working array.                                    *
*    (8) z : 1-dim. working array.                                     *
*    (9) ier : error code.                                             *
*  copy right   t. oguni   july 30 1989   version 1.0                  *
************************************************************************
	integer n, nud, n1, np, nq, ier
	complex*16 a(n1,n), b(n), dr(n), z(n)
	real*8 eps
	complex*16 xx, s, sum, au, t
	real*8 eps1
	integer i ,m, j, k1, mj, i1, k, j1
c  check the input data
	ier = 0
	eps1 = 1.0d-14
	m = nud + 1
	if ((n .le. 0) .or. (nud .le. 0 ) .or. (n1 .lt. m)) then
	 ier = 2
	 write(*,*) '(subr. lisb) invalid argument. ', n, nud, n1
	 return
	endif
	if (eps .le. 0.0) eps = eps1
c  modified cholesky decomposition
	j = 1
	if (cdabs(a(m,1)) .le. eps) then
	 ier = 1
	 write(*,*) '(subr. lisb) singular at step # ', j
	 return
	endif
	dr(1) = dcmplx(1.0d0) / a(m,1)
	xx = a(m-1,2)
	a(m-1,2) = a(m-1,2) * dr(1)
	s = a(m,2) - xx * a(m-1,2)
	j = 2
	if (cdabs(s) .le. eps) then
	 ier = 1
	 write(*,*) '(subr. lisb) singular at step # ', j
	 return
	endif
	dr(2) = dcmplx(1.0d0) / s
	if (m .lt. 3) then
	 do 5 j=3,n
	  xx = a(1,j)
	  a(1,j) = xx * dr(j-1)
	  s = a(2,j) - xx * a(1,j)
	  if (cdabs(s) .le. eps) then
	   ier = 1
	   write(*,*) ' (subr. lisb) singular at step # ', j
	   return
	  endif
	  dr(j) = dcmplx(1.0d0) / s
    5	 continue
	else
	 do 30 j=3,n
	  k1 = 1
	  if (j .ge. m) k1 = j - m + 1
	  mj = m - j
	  do 20 i=k1+1,j-1
	   sum = dcmplx(0.0d0)
	   do 10 k=k1,i-1
   10	    sum = sum + a(m-i+k,i) * a(mj+k,j)
	   a(mj+i,j) = a(mj+i,j) - sum
   20	  continue
	  sum = dcmplx(0.0d0)
	  do 25 i=k1,j-1
	   xx = a(mj+i,j)
	   au = xx * dr(i)
	   sum = sum + xx *au
	   a(mj+i,j) = au
   25	  continue
	  t = a(m,j) - sum
	  if (cdabs(t) .le. eps) then
	   ier = 1
	   write(*,*) ' (subr. lisb) singular at step # ', j
	   return
	  endif
	  dr(j) = dcmplx(1.0d0) / t
   30	 continue
	endif
c subtitution
	entry dcsbsub(a, n, nud, n1, np, nq, b, eps, dr, z, ier)
c  forward substitution
	m = nud + 1
	if (m .lt. 3) then
	  z(np) = b(np)
	  do 40 j=np+1,n
   40	  z(j) = b(j) - a(1,j) * z(j-1)
	  do 45 j=1,np-1
   45	  z(j) = dcmplx(0.d0)
	  do 50 j=np,n
   50	  z(j) = z(j) * dr(j)
	  b(n) = z(n)
c	  do 60 j=1,n-1
	  do 60 j=1,nq-1
   60	  b(n-j) = z(n-j) - a(1,n-j+1) * b(n-j+1)
	else
	  z(np) = b(np)
	  z(np+1) = b(np+1) - a(m-1,np+1) * z(np)
	  do 80 j=np+2,n
	    if (j .gt. np-1+m) then
	      i1 = 1
	    else
	      i1 = np-1+m - j + 1
	    endif
	    sum = dcmplx(0.0d0)
	    do 70 k=i1,m-1
   70       sum = sum + a(k,j) * z(j-m+k)
   80	  z(j) = b(j) - sum
	  do 85 j=1,np-1
   85	  z(j) = dcmplx(0.d0)
	  do 90 j=np,n
   90	  z(j) = z(j) * dr(j)
c
	  b(n) = z(n)
	  b(n-1) = z(n-1) - a(m-1,n) * z(n)
	  do 110 j=3,nq
	    j1 = n - j + 1
	    i1 = 1
	    if (j .lt. m) i1 = m - j + 1
	    sum = dcmplx(0.0d0)
	    do 100 k=i1,m-1
  100	    sum = sum + a(k,m-k+j1) * b(m-k+j1)
	    b(j1) = z(j1) - sum
  110	  continue
	endif
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine init_complex_array
     &	             ( n_station,station_displacement )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c initializing the acculuated displacement at the station.
c    required subroutines: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
	integer n_station
	complex*16 station_displacement(*)
c other variables
	integer i_station
c
	do 100 i_station=1,n_station
	  station_displacement(i_station) = dcmplx(0.d0)
  100	continue
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine comp_displacement_station
     &	    ( maxngrid_r,maxlmax,whole_vector_tor,whole_vector_sph,
     &	      l,n_station,station_theta,station_phi,
     &	      idim_station_sph,idim_station_tor,
     &	      vecsph_sph1,vecsph_sph2,vecsph_tor,
     &	      station_displacement )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c accumulating the displacement at the station.
c    required subroutines: error_handling,comp_vecsph_tor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
	integer maxngrid_r,maxlmax,l,n_station
	integer idim_station_sph,idim_station_tor
	real*8 station_theta(*),station_phi(*)
	complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
	complex*16 whole_vector_tor(maxngrid_r,-2:2)
	complex*16 vecsph_sph1(3,0:maxlmax,-2:2,*)
	complex*16 vecsph_sph2(3,0:maxlmax,-2:2,*)
	complex*16 vecsph_tor(3,0:maxlmax,-2:2,*)
	complex*16 station_displacement(3,*)
c other variables
	integer m,i_station,icomp
	real*8 lsq
c
	if ( l.le.0 ) call error_handling(55)
	lsq = dsqrt( dble(l) * dble(l+1) )
	do 200 i_station=1,n_station
c ---- computing the value of the trial functions at the station
	  do 120 m=max0(-l,-2),min0(l,2)
c -------- horizontal dependent part
c	    call comp_vecsph(l,m,
c     &	                     station_theta(i_station),
c     &	                     station_phi(i_station),
c     &	                     vecsph_sph1,vecsph_sph2,vecsph_tor)
c ---- computing the displacement at the station
	    do 110 icomp=1,3
	      station_displacement(icomp,i_station)
     &	      = station_displacement(icomp,i_station)
     &	        + whole_vector_sph(idim_station_sph,m)
     &	          * vecsph_sph1(icomp,l,m,i_station)
     &	        + whole_vector_sph(idim_station_sph+1,m)
     &	          * vecsph_sph2(icomp,l,m,i_station) / dcmplx( lsq )
     &	        + whole_vector_tor(idim_station_tor,m)
     &	          * vecsph_tor(icomp,l,m,i_station) / dcmplx( lsq )
  110	    continue
  120	  continue
  200	continue
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine comp_displacement_station0
     &	    ( maxngrid_r,maxlmax,whole_vector_tor,whole_vector_sph,
     &	      l,n_station,station_theta,station_phi,
     &	      idim_station_sph,idim_station_tor,
     &	      vecsph_sph1,vecsph_sph2,vecsph_tor,
     &	      station_displacement )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c accumulating the displacement at the station.
c    required subroutines: error_handling,comp_vecsph_tor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
	integer maxngrid_r,maxlmax,l,n_station
	integer idim_station_sph,idim_station_tor
	real*8 station_theta(*),station_phi(*)
	complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
	complex*16 whole_vector_tor(maxngrid_r,-2:2)
	complex*16 vecsph_sph1(3,0:maxlmax,-2:2,*)
	complex*16 vecsph_sph2(3,0:maxlmax,-2:2,*)
	complex*16 vecsph_tor(3,0:maxlmax,-2:2,*)
	complex*16 station_displacement(3,*)
c other variables
	integer m,i_station,icomp
	real*8 lsq
c
	if ( l.ne.0 ) call error_handling(56)
	lsq = dsqrt( dble(l) * dble(l+1) )
	do 200 i_station=1,n_station
c ---- computing the value of the trial functions at the station
	  do 120 m=max0(-l,-2),min0(l,2)
c -------- horizontal dependent part
c	    call comp_vecsph(l,m,
c     &	                     station_theta(i_station),
c     &	                     station_phi(i_station),
c     &	                     vecsph_sph1,vecsph_sph2,vecsph_tor)
c ---- computing the displacement at the station
	    do 110 icomp=1,3
	      station_displacement(icomp,i_station)
     &	      = station_displacement(icomp,i_station)
     &	        + whole_vector_sph(idim_station_sph,m)
     &	          * vecsph_sph1(icomp,l,m,i_station)
  110	    continue
  120	  continue
  200	continue
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       
        subroutine comp_disp_strain_all(station_disp_solid,
     &            station_disp_fluid,
     &            epsilon11,epsilon22,epsilon33,epsilon12,
     &            epsilon13,epsilon23,pressure,
     &            lmax,coef_c,coef_dcdr,maxngrid_r,max_ndep,
     &            r_sta_solid,r_sta_fluid,ir_dep_solid,ir_dep_fluid,
     &            ir_for_stress,ir_for_pressure,
     &            ir_top_stress,ir_bot_stress,
     &            ir_top_pressure,ir_bot_pressure,
     &            idim_ir_sph,idim_ir_tor,
     &            whole_vector_tor,whole_vector_sph,
     &            nl_check_amp,i_significance,l,ir_CMB,
     &            ir_ICB,grid_rho,n_station_solid,n_station_fluid,
     &            ndist_ela,ndist_acou,dist_elastic,dist_acoustic,rank,
     &            omega,grid_r,ngrid_r,depth_solid,depth_fluid,
     &            ndep_solid,ndep_fluid,top_fluid,bot_fluid,
     &            source_type,ifreq,pi)

c input/output
        implicit none
        integer ifreq,lmax
        complex*16 station_disp_solid(3,*)
        complex*16 station_disp_fluid(3,*)
        complex*16 epsilon11(*),epsilon22(*),
     &            epsilon33(*),epsilon12(*),
     &            epsilon13(*),epsilon23(*)
        complex*16 pressure(*)
        complex*16 coef_c(0:lmax,-2:2,2,3),coef_dcdr(0:lmax,-2:2,2,3)
        integer    maxngrid_r,max_ndep,l,n_station_solid
        integer  nl_check_amp,i_significance
        integer  n_station_fluid,ngrid_r
        real*8   r_sta_solid(max_ndep),r_sta_fluid(max_ndep)
        integer    ir_dep_solid(max_ndep),ir_dep_fluid(max_ndep)
        integer    ir_for_stress(3,max_ndep)
        integer    ir_for_pressure(3,max_ndep)
        integer  ir_top_stress(3),ir_bot_stress(3)
        integer  ir_top_pressure(3),ir_bot_pressure(3)

        integer idim_ir_sph(maxngrid_r),idim_ir_tor(maxngrid_r)
        complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
        complex*16 whole_vector_tor(maxngrid_r,-2:2)
        complex*16 omega
        integer ir_CMB,ir_ICB
        real*8 grid_rho(2,maxngrid_r)
        real*8 depth_solid(max_ndep),depth_fluid(max_ndep)
c        real*8  min_dep_solid,min_dep_fluid,ddep
        integer ndep_solid,ndep_fluid
        integer ndist_ela,ndist_acou
        real*8  dist_acoustic(*),dist_elastic(*)
        integer rank
        real*8  grid_r(maxngrid_r)
        integer source_type
        integer top_fluid,bot_fluid
        real*8  pi

c other variables
        integer idep,idep_new,itheta,ntheta
        integer istation_solid,istation_fluid,ir
        integer ir_up_stress,ir_mid_stress,ir_down_stress
        integer ir_up_pressure,ir_mid_pressure,ir_down_pressure
        integer m,max_m,min_m
        real*8  theta,phi,depth,r_station
        complex*16 ur,utheta,uphi
        complex*16 epsilon11_lm,epsilon22_lm,epsilon33_lm
        complex*16 epsilon12_lm,epsilon13_lm,epsilon23_lm
        complex*16 pressure_lm
        complex*16 coef_c_lm(3),coef_dcdr_lm(3)
        logical fluid,copy_top_coef,copy_bot_coef
        integer m_cycles,nl_for_average
        real*8  dl_onecycle

        if(source_type.eq.1) then
          max_m=2
          min_m=-2
        else if(source_type.eq.2) then
          max_m=1
          min_m=-1
        end if

            
        if(ndep_solid.gt.0) then
           ntheta=ndist_ela
        else
           ntheta=ndist_acou
        end if
        do 100 itheta=1,ntheta 
         if(ndep_solid.gt.1) then
           theta=dist_elastic(itheta)
         else
           theta=dist_acoustic(itheta)
         end if

         phi=0.0
!        In order to remove the undulation for grids around the source depth, we take the average values within M periods.
!        The undulation period of sum(c_lm*Ylm) (when l>>kr and k is the wavenumber) equals dl_onecycle=2*PI/theta. Thus, We take
!        nl_for_average as M*dl_period, where M is the largest integer 
!        making M*dl_onecycle smaller than nl_check_amp.

         dl_onecycle=2*pi/theta
         m_cycles=int((nl_check_amp/dl_onecycle))
         nl_for_average=int(m_cycles*dl_onecycle)
         if(nl_for_average.lt.0.or.nl_for_average.gt.nl_check_amp) then
           write(6,*) 'Error, nl_for_average <0 or >nl_check_amp'
           stop
         end if
         
         if(nl_for_average.lt.100) then
          write(6,*) "WARNING:nl_for_average<100!!!!!!! 
     &         Truncation error during removing undulation 
     &          might be as large as 1/100. Possible reasons might be
     &          nl_check_amp < 100 or distance is not teleseismic
     &          (<3.6deg). Try nl_check_amp=>150 or get rid of station
     &          with distance smaller than 3.6deg."
         end if

!We want to calculate some variables, such as Ylm, only once when idep=1(see sub-subfuction comp_disp_strain), 
!so solid and fluid media are packaged in one loop instead of two loops.
         do 200 idep=1,ndep_solid+ndep_fluid
!fluid media
           if(idep>ndep_solid) then
              fluid=.TRUE.
              idep_new=idep-ndep_solid
c              depth=min_dep_fluid+(idep_new-1)*ddep
              depth=depth_fluid(idep_new)
              r_station=r_sta_fluid(idep_new)
              ir=ir_dep_fluid(idep_new)
              ir_up_pressure=ir_for_pressure(3,idep_new)
              ir_mid_pressure=ir_for_pressure(2,idep_new)
              ir_down_pressure=ir_for_pressure(1,idep_new)
              istation_fluid=(idep_new-1)*ntheta+itheta
              istation_solid=1
!solid
           else
              fluid=.FALSE.
              idep_new=idep
c              depth=min_dep_solid+(idep_new-1)*ddep
              depth=depth_solid(idep_new)
              r_station=r_sta_solid(idep_new)
              ir=ir_dep_solid(idep_new)
              ir_up_stress=ir_for_stress(3,idep_new)
              ir_mid_stress=ir_for_stress(2,idep_new)
              ir_down_stress=ir_for_stress(1,idep_new)
              istation_solid=(idep_new-1)*ntheta+itheta
              istation_fluid=1
           end if
 
           copy_top_coef=.false.
           copy_bot_coef=.false.
           if(top_fluid.eq.1.and.fluid.and.idep.eq.ndep_solid+1) 
     &          copy_top_coef=.true.
           if(top_fluid.eq.0.and..not.fluid.and.idep.eq.1) 
     &          copy_top_coef=.true.
           if(bot_fluid.eq.1.and.fluid.and.
     &        idep.eq.ndep_solid+ndep_fluid) 
     &          copy_bot_coef=.true.
           if(bot_fluid.eq.0.and..not.fluid.and.idep.eq.ndep_solid) 
     &          copy_bot_coef=.true.

           do 120 m=max0(-l,min_m),min0(l,max_m)
               
              call comp_disp_strain(omega,l,m,theta,phi,ir,
     &           ngrid_r,grid_r,ifreq,ir_up_stress,
     &           ir_mid_stress,ir_down_stress,
     &           ir_up_pressure,ir_mid_pressure,ir_down_pressure,
     &           ir_top_stress,ir_bot_stress,
     &           ir_top_pressure,ir_bot_pressure,
     &           ir_CMB,ir_ICB,grid_rho,whole_vector_sph,
     &           whole_vector_tor,maxngrid_r,
     &           fluid,idim_ir_sph,idim_ir_tor,idep,idep_new,
     &           ndep_solid,ndep_fluid,r_station,
     &           ur,utheta,uphi,pressure_lm,epsilon11_lm,epsilon22_lm,
     &           epsilon33_lm,epsilon12_lm,epsilon13_lm,epsilon23_lm,
     &           coef_c_lm,coef_dcdr_lm,top_fluid,bot_fluid)

              call sum_disp_strain(station_disp_solid(1,istation_solid),
     &           station_disp_fluid(1,istation_fluid),fluid,ur,utheta,
     &           uphi,istation_solid,istation_fluid,epsilon11,
     &           epsilon22,epsilon33,epsilon12,epsilon13,epsilon23,
     &           epsilon11_lm,epsilon22_lm,epsilon33_lm,epsilon12_lm,
     &           epsilon13_lm,epsilon23_lm,pressure,pressure_lm,
     &           copy_top_coef,copy_bot_coef,coef_c_lm,coef_dcdr_lm,
     &           l,m,nl_check_amp,i_significance,lmax,coef_c,coef_dcdr,
     &           nl_for_average)
!for convegence test
!           if(ifreq.eq.50.and.m.eq.0.and.idep_new.eq.1.and.itheta.eq.1) print
!          if(ifreq.eq.50.and.m.eq.0)
!     &     print *,'sum=',
!     &          station_disp_solid(1,istation_solid),l,idep_new

  120      continue
  200    continue

  100  continue
        
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine comp_disp_strain_all0(station_disp_solid,
     &            station_disp_fluid,
     &            epsilon11,epsilon22,epsilon33,epsilon12,
     &            epsilon13,epsilon23,pressure,
     &            nl_check_amp,i_significance,
     &            lmax,coef_c,coef_dcdr,maxngrid_r,max_ndep,
     &            r_sta_solid,r_sta_fluid,
     &            ir_dep_solid,ir_dep_fluid,idim_ir_sph0,
     &            ir_for_stress,ir_for_pressure,
     &            ir_top_stress,ir_bot_stress,
     &            ir_top_pressure,ir_bot_pressure,
     &            whole_vector_sph,l,ir_CMB,ir_ICB,grid_rho,
     &            n_station_solid,n_station_fluid,ndist_ela,
     &            ndist_acou,dist_elastic,dist_acoustic,
     &            omega,grid_r,ngrid_r,
     &            depth_solid,depth_fluid,
     &            ndep_solid,ndep_fluid,top_fluid,bot_fluid,
     &            source_type,ifreq)
        implicit none
c input/output
        complex*16 station_disp_solid(3,*),station_disp_fluid(3,*)
        complex*16 epsilon11(*),epsilon22(*),
     &             epsilon33(*),epsilon12(*),
     &             epsilon13(*),epsilon23(*)
        complex*16 pressure(*)
        complex*16 coef_c(0:lmax,-2:2,2,3),coef_dcdr(0:lmax,-2:2,2,3)
        integer    maxngrid_r,max_ndep,l,lmax,ifreq
        integer  nl_check_amp,i_significance
        integer  n_station_solid,n_station_fluid,ngrid_r
        integer  ir_dep_solid(max_ndep),ir_dep_fluid(max_ndep)
        integer  ir_for_stress(3,max_ndep)
        integer  ir_for_pressure(3,max_ndep)
        integer  ir_top_stress(3),ir_bot_stress(3)
        integer  ir_top_pressure(3),ir_bot_pressure(3)
        real*8  r_sta_solid(max_ndep),r_sta_fluid(max_ndep)
        integer idim_ir_sph0(maxngrid_r)
        complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
        complex*16 omega
        integer ir_CMB,ir_ICB
        real*8 grid_rho(2,maxngrid_r)
        real*8  depth_solid(max_ndep),depth_fluid(max_ndep)
c        real*8  min_dep_solid,min_dep_fluid,ddep
        integer ndep_solid,ndep_fluid
        integer ndist_acou,ndist_ela
        real*8  dist_acoustic(*),dist_elastic(*)
        integer source_type
        integer top_fluid,bot_fluid
        real*8  grid_r(maxngrid_r)

c other variables
        integer ir,idep,itheta,ntheta,idep_new
        integer istation_solid,istation_fluid
        integer ir_up_stress,ir_mid_stress,ir_down_stress
        integer ir_up_pressure,ir_mid_pressure,ir_down_pressure
        real*8  theta,phi,depth,r_station
        complex*16 ur,utheta,uphi
        complex*16 epsilon11_lm,epsilon22_lm,epsilon33_lm
        complex*16 epsilon12_lm,epsilon13_lm,epsilon23_lm
        complex*16 pressure_lm
        complex*16 coef_c_lm(3),coef_dcdr_lm(3)
        logical fluid,copy_top_coef,copy_bot_coef
        real*8 degtorad
        parameter ( degtorad=3.1415926535897932d0/180.d0)



        do 100 idep=1,ndep_solid+ndep_fluid
          !solid media
          if(idep>ndep_solid) then
              fluid=.TRUE.
              idep_new=idep-ndep_solid
c              depth=min_dep_fluid+(idep_new-1)*ddep
              depth=depth_fluid(idep_new)
              r_station=r_sta_fluid(idep_new)
              ir=ir_dep_fluid(idep_new)
              ir_up_pressure=ir_for_pressure(3,idep_new)
              ir_mid_pressure=ir_for_pressure(2,idep_new)
              ir_down_pressure=ir_for_pressure(1,idep_new)
              ntheta=ndist_acou
              istation_fluid=(idep_new-1)*ntheta
              istation_solid=1
!fluid
          else
              fluid=.FALSE.
              idep_new=idep
c              depth=min_dep_solid+(idep_new-1)*ddep
              depth=depth_solid(idep_new)
              r_station=r_sta_solid(idep_new)
              ir=ir_dep_solid(idep_new)
              ir_up_stress=ir_for_stress(3,idep_new)
              ir_mid_stress=ir_for_stress(2,idep_new)
              ir_down_stress=ir_for_stress(1,idep_new)
              ntheta=ndist_ela
              istation_solid=(idep_new-1)*ntheta
              istation_fluid=1
          end if

           copy_top_coef=.false.
           copy_bot_coef=.false.
           if(top_fluid.eq.1.and.fluid.and.idep.eq.ndep_solid+1)
     &          copy_top_coef=.true.
           if(top_fluid.eq.0.and..not.fluid.and.idep.eq.1)
     &          copy_top_coef=.true.
           if(bot_fluid.eq.1.and.fluid.and.
     &        idep.eq.ndep_solid+ndep_fluid)
     &          copy_bot_coef=.true.
           if(bot_fluid.eq.0.and..not.fluid.and.idep.eq.ndep_solid)
     &          copy_bot_coef=.true.

         do 200 itheta=1,ntheta
               if(fluid) then
                  istation_fluid=istation_fluid+1
                  theta=dist_acoustic(itheta)
               else
                  istation_solid=istation_solid+1
                  theta=dist_elastic(itheta)
               end if
               phi=0.0

               call comp_disp_strain0(omega,l,0,theta,phi,ir,
     &           ngrid_r,grid_r,ifreq,ir_up_stress,
     &           ir_mid_stress,ir_down_stress,
     &           ir_up_pressure,ir_mid_pressure,ir_down_pressure,
     &           ir_top_stress,ir_bot_stress,
     &           ir_top_pressure,ir_bot_pressure,
     &           ir_CMB,ir_ICB,whole_vector_sph,maxngrid_r,
     &           fluid,idim_ir_sph0,idep,idep_new,ndep_solid,
     &           ndep_fluid,grid_rho,r_station,
     &           ur,utheta,uphi,pressure_lm,epsilon11_lm,epsilon22_lm,
     &           epsilon33_lm,epsilon12_lm,epsilon13_lm,epsilon23_lm,
     &           coef_c_lm,coef_dcdr_lm,top_fluid,bot_fluid)

               call sum_disp_strain(
     &           station_disp_solid(1,istation_solid),
     &           station_disp_fluid(1,istation_fluid),fluid,ur,utheta,
     &           uphi,istation_solid,istation_fluid,epsilon11,
     &           epsilon22,epsilon33,epsilon12,epsilon13,epsilon23,
     &           epsilon11_lm,epsilon22_lm,epsilon33_lm,epsilon12_lm,
     &           epsilon13_lm,epsilon23_lm,pressure,pressure_lm,
     &           copy_top_coef,copy_bot_coef,coef_c_lm,coef_dcdr_lm,
     &           l,0,nl_check_amp,i_significance,lmax,coef_c,coef_dcdr,
     &           1)

  200    continue
  100   continue

        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine comp_disp_strain(omega,l,m,theta,phi,
     &           ir,ngrid_r,grid_r,ifreq,ir_up_stress,
     &           ir_mid_stress,ir_down_stress,
     &           ir_up_pressure,ir_mid_pressure,ir_down_pressure,
     &           ir_top_stress,ir_bot_stress,
     &           ir_top_pressure,ir_bot_pressure,
     &           ir_CMB,ir_ICB,grid_rho,whole_vector_sph,
     &           whole_vector_tor,maxngrid_r,
     &           fluid,idim_ir_sph,idim_ir_tor,idep,idep_new,
     &           ndep_solid,ndep_fluid,r_station,
     &           ur,utheta,uphi,pressure_lm,epsilon11_lm,epsilon22_lm,
     &           epsilon33_lm,epsilon12_lm,epsilon13_lm,epsilon23_lm,
     &           coef_c_lm,coef_dcdr_lm,top_fluid,bot_fluid)

        implicit none
c variables for input/output
        integer ifreq
        integer maxngrid_r,l,m,ir,ngrid_r,idep,idep_new,ir_CMB,ir_ICB
        integer ndep_solid,ndep_fluid,top_fluid,bot_fluid
        integer ir_up_stress,ir_mid_stress,ir_down_stress
        integer ir_up_pressure,ir_mid_pressure,ir_down_pressure
        integer  ir_top_stress(3),ir_bot_stress(3)
        integer  ir_top_pressure(3),ir_bot_pressure(3)
        real*8 theta,phi,r_station
        real*8 grid_r(*)
        complex*16 omega
        complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
        complex*16 whole_vector_tor(maxngrid_r,-2:2)
        real*8 grid_rho(2,maxngrid_r)
        logical fluid
        integer idim_ir_sph(maxngrid_r),idim_ir_tor(maxngrid_r)

        complex*16 ur,utheta,uphi
        complex*16 epsilon11_lm,epsilon22_lm,epsilon33_lm
        complex*16 epsilon12_lm,epsilon13_lm,epsilon23_lm
        complex*16 pressure_lm
        complex*16 coef_c_lm(3),coef_dcdr_lm(3)

c other variables 
        integer m0,i
        real*8 factor(-2:2),x,y
        complex*16 c1,dc1,c2,dc2,c3,dc3
        real*8 dr
        real*8 lsq
        real*8 plm(-2:2),plm1(-2:2),plm2(-2:2)
        real*8 plgndr
        real*8 dpdtheta(-2:2),d2pd2theta(-2:2)
        real*8 dplm1dtheta(-2:2)
        real*8 rho,rho_up,rho_down
        real*8 r_input(3)
        complex*16 c_input(9),c_out_interpo(3)
        complex*16 dcdr_out_interpo(3)
        complex*16 dphi(-2:2),Ylm(-2:2),dYdphi(-2:2)
        complex*16 dYdtheta(-2:2),d2Yd2theta(-2:2)
        complex*16 dc1dr,dc2dr,dc3dr
        complex*16 du1dr,du1dtheta,du1dphi
        complex*16 du2dr,du2dtheta,du2dphi
        complex*16 du3dr,du3dtheta,du3dphi
        complex*16 dQdr,dQdtheta,dQdphi
        complex*16 expimp(-2:2)
        integer idim_sph,idim_tor

        integer idim_up_pressure,idim_down_pressure
        integer idim_mid_pressure
        integer idim_up_S_stress,idim_up_T_stress
        integer idim_mid_S_stress,idim_mid_T_stress
        integer idim_down_S_stress,idim_down_T_stress
        integer idim1

        save factor,x,y,lsq,plm,plm1,plm2,dpdtheta,d2pd2theta
        save dplm1dtheta,dphi,Ylm,dYdphi,dYdtheta,d2Yd2theta,expimp
c constant 
        integer ncomp,icomp
        real*8 pi
        parameter ( ncomp= 3)
        parameter ( pi=3.1415926535897932d0)


        
        m0=iabs(m)
        if ( (l.lt.0).or.(m0.gt.l) ) call error_handling(41)
c **********************************************************************
c computing the normalization factor (including the sign)
c **********************************************************************
c
        if(idep.eq.1) then
         factor(m) = 1.d0
         do 100 i=l-m0+1,l+m0
           factor(m) = factor(m) * dble(i)
  100    continue
         factor(m) = dsqrt( dble(2*l+1)/(4.d0*pi) / factor(m) )
         if ( ( m0.ne.m ).and.( mod(m0,2).eq.1) ) factor(m) = -factor(m)
         lsq=dsqrt(dble(l)*dble(l+1))
c **********************************************************************
c computing the displacement and its derivation
c **********************************************************************
         x=dcos(theta)
         y=dsin(theta)
         expimp(m)=cdexp( dcmplx( 0.d0, dble(m)*phi ) )
         plm(m)=plgndr(l,m0,x)
         if(l.eq.1) then
            if(m0.eq.1) then
              dpdtheta(m)=-x
              d2pd2theta(m)=y
            else if(m0.eq.0) then
              dpdtheta(m)=-y
              d2pd2theta(m)=-x
            end if
         else if(l.eq.2.and.m0.eq.2) then
            dpdtheta(m)=6.d0*x*y
            d2pd2theta(m)=6.d0*(x*x-y*y)
         else 
            plm1(m)=plgndr(l,m0+1,x)
            dpdtheta(m)=(dble(m0)*x/y*plm(m)+plm1(m))
            if(m0+2.le.l) then
              plm2(m)=plgndr(l,m0+2,x)
              dplm1dtheta(m)=(dble(m0+1)*x/y*plm1(m)+plm2(m))
            else
              dplm1dtheta(m)=-dble(l+m0+1)*dble(l-m0)*plm(m)
     &                     -dble(m0+1)*x/y*plm1(m)
            end if
            d2pd2theta(m)=dble(m0)*x/y*(-plm(m)/(y*x)
     &                    + dpdtheta(m))+dplm1dtheta(m)
         end if
         dphi(m)=dcmplx(0.d0,dble(m))
         Ylm(m)=factor(m)*expimp(m)*plm(m)
         dYdtheta(m)=factor(m)*expimp(m)*dpdtheta(m)
         d2Yd2theta(m)=factor(m)*expimp(m)*d2pd2theta(m)
         dYdphi(m)=Ylm(m)*dphi(m)
        end if
c***********************************************************************
c computing the displacemetn
c***********************************************************************

        c1=dcmplx(0.d0,0.d0)
        c2=dcmplx(0.d0,0.d0)
        c3=dcmplx(0.d0,0.d0)
        dc1dr=dcmplx(0.d0,0.d0)
        dc2dr=dcmplx(0.d0,0.d0)
        dc3dr=dcmplx(0.d0,0.d0)

c***********************fluid media***********************************
        if(fluid) then
! need to test rho
         if(ir.eq.ir_CMB) then
             idim_sph=idim_ir_sph(ir-1)+1
             rho=grid_rho(2,ir-1)
         else if(ir.eq.ir_CMB-1) then
             idim_sph=idim_ir_sph(ir)
             rho=grid_rho(1,ir)
         else if(ir.lt.ir_CMB-1.and.ir.gt.ir_ICB) then
!Central differencing scheme
             idim_sph=idim_ir_sph(ir)+1
             rho=grid_rho(1,ir)
         else if(ir.eq.ir_ICB) then
             idim_sph=idim_ir_sph(ir)+1
             rho=grid_rho(1,ir)
         else  
c             print *,'ir_ICB',ir_ICB,ir,grid_r(ir)
             stop 'Error in finding idim_sph in the outer core'
         end if
!idim for ir_up_pressure
         if(ir_up_pressure.eq.ir_CMB) then
             idim_up_pressure=idim_ir_sph(ir_up_pressure-1)+1
         else if(ir_up_pressure.eq.ir_CMB-1) then
             idim_up_pressure=idim_ir_sph(ir_up_pressure)
         else if(ir_up_pressure.lt.ir_CMB-1.and.
     &           ir_up_pressure.gt.ir_ICB) then
             idim_up_pressure=idim_ir_sph(ir_up_pressure)+1
         else if(ir_up_pressure.eq.ir_ICB) then
             idim_up_pressure=idim_ir_sph(ir_up_pressure)+1
         else
            stop 'Error in finding idim_sph_up in the outer core'
         end if

         if(ir_mid_pressure.eq.ir_CMB) then
             idim_mid_pressure=idim_ir_sph(ir_mid_pressure-1)+1
         else if(ir_mid_pressure.eq.ir_CMB-1) then
             idim_mid_pressure=idim_ir_sph(ir_mid_pressure)
         else if(ir_mid_pressure.lt.ir_CMB-1.and.
     &           ir_mid_pressure.gt.ir_ICB) then
             idim_mid_pressure=idim_ir_sph(ir_mid_pressure)+1
         else if(ir_mid_pressure.eq.ir_ICB) then
             idim_mid_pressure=idim_ir_sph(ir_mid_pressure)+1
         else
            stop 'Error in finding idim_sph_mid in the outer core'
         end if

         !idim for ir_down_pressure
         if(ir_down_pressure.eq.ir_CMB) then
             idim_down_pressure=idim_ir_sph(ir_down_pressure-1)+1
         else if(ir_down_pressure.eq.ir_CMB-1) then
             idim_down_pressure=idim_ir_sph(ir_down_pressure)
         else if(ir_down_pressure.lt.ir_CMB-1.and.
     &           ir_down_pressure.gt.ir_ICB) then
             idim_down_pressure=idim_ir_sph(ir_down_pressure)+1
         else if(ir_down_pressure.eq.ir_ICB) then
             idim_down_pressure=idim_ir_sph(ir_down_pressure)+1
         else
           stop 'Error in finding idim_sph_down in the outer core'
         end if

         c1=whole_vector_sph(idim_sph,m)
         pressure_lm=omega*c1*Ylm(m)

         dc1=whole_vector_sph(idim_up_pressure,m)-
     &      whole_vector_sph(idim_down_pressure,m)
         dr=grid_r(ir_up_pressure)-grid_r(ir_down_pressure)
         dc1dr=dc1/dr
!         if(ir_up_pressure.eq.ir_CMB) then
!           rho_up=grid_rho(2,ir_up_pressure-1)
!         else
!           rho_up=grid_rho(1,ir_up_pressure)
!         end if

!         if(ir_mid_pressure.eq.ir_CMB) then
!           rho_mid=grid_rho(2,ir_mid_pressure-1)
!         else
!           rho_mid=grid_rho(1,ir_mid_pressure)
!         end if


!         if(ir_down_pressure.eq.ir_CMB) then
!           rho_down=grid_rho(2,ir_down_pressure-1)
!         else
!           rho_down=grid_rho(1,ir_down_pressure)
!         end if


         dQdr=dc1*Ylm(m)/dr
         dQdtheta=c1*dYdtheta(m)
         dQdphi=c1*dYdphi(m)
         ur=-dQdr/(rho*omega)
         utheta=-dQdtheta/(rho*omega*grid_r(ir))
c         uphi=-dQdphi/(rho*omega*grid_r(ir)*sin(theta))
         uphi=dcmplx(0.d0,0.d0)
cc         uphi=omega*c1*Ylm(m)

c ***********************solid media***********************************
        else
!for quadra inperpolation
         if(ir.eq.ir_ICB) then
             idim_sph=idim_ir_sph(ir-1)
         else if(ir.eq.ir_ICB-1) then
             idim_sph=idim_ir_sph(ir)-2
ccccc    idim_ir_sph(ir-1) is equal to idim_ir_sph(ir)-4
         else if(ir.gt.3) then
             idim_sph=idim_ir_sph(ir)
             if(ir.ge.ir_CMB)
     &         idim_tor=idim_ir_tor(ir)
         else 
c             print *,'check_ir_idep',ir,idep_new
             stop 'Error in finding targeted idim_ir'
         end if

         if(ir_up_stress.eq.ir_ICB) then
            idim_up_S_stress=idim_ir_sph(ir_up_stress-1)
         else if(ir_up_stress.eq.ir_ICB-1) then
            idim_up_S_stress=idim_ir_sph(ir_up_stress)-2
         else if(ir_up_stress.gt.3) then
            idim_up_S_stress=idim_ir_sph(ir_up_stress)
            if(ir_up_stress.ge.ir_CMB) 
     &        idim_up_T_stress=idim_ir_tor(ir_up_stress)
         else
            stop 'Error in finding targeted idim_up_stress'
         end if

         if(ir_mid_stress.eq.ir_ICB) then
            idim_mid_S_stress=idim_ir_sph(ir_mid_stress-1)
         else if(ir_mid_stress.eq.ir_ICB-1) then
            idim_mid_S_stress=idim_ir_sph(ir_mid_stress)-2
         else if(ir_mid_stress.gt.3) then
            idim_mid_S_stress=idim_ir_sph(ir_mid_stress)
            if(ir_mid_stress.ge.ir_CMB)
     &        idim_mid_T_stress=idim_ir_tor(ir_mid_stress)
         else
            stop 'Error in finding targeted idim_up_stress'
         end if

            
         if(ir_down_stress.eq.ir_ICB) then
            idim_down_S_stress=idim_ir_sph(ir_down_stress-1)
         else if(ir_down_stress.eq.ir_ICB-1) then
            idim_down_S_stress=idim_ir_sph(ir_down_stress)-2
         else if(ir_down_stress.gt.3) then
            idim_down_S_stress=idim_ir_sph(ir_down_stress)
            if(ir_down_stress.ge.ir_CMB)
     &        idim_down_T_stress=idim_ir_tor(ir_down_stress)
         else
c            print *,'ir_down_stress',ir_down_stress,idep_new
            stop 'Error in finding targeted idim_down_stress'
         end if


c***********************************************************************
c computing displacement and strain
c***********************************************************************
!WENBO to fix

         c_input=dcmplx(0.d0,0.d0)
         do i=1,3
          do icomp=1,ncomp
!The most up grid point between ir,ir_up,ir_down
            if(i.eq.1) then
                 r_input(i)=grid_r(ir_up_stress)
!c1
                 if(icomp.eq.1) then
                   c_input( ncomp * (i-1) + icomp)=
     &               whole_vector_sph(idim_up_S_stress,m)
!c2
                  else if(icomp.eq.2) then
                   c_input( ncomp * (i-1) + icomp)=
     &               whole_vector_sph(idim_up_S_stress+1,m)
!c3
                 else if(icomp.eq.3.and.ir.ge.ir_CMB) then
                   c_input( ncomp * (i-1) + icomp)=
     &               whole_vector_tor(idim_up_T_stress,m)
                 end if
!the middle grid point
            else if(i.eq.2) then
               r_input(i)=grid_r(ir_mid_stress)
               if(icomp.eq.1) then
                 c_input( ncomp * (i-1) + icomp)=
     &             whole_vector_sph(idim_mid_S_stress,m)
               else if(icomp.eq.2) then
                 c_input( ncomp * (i-1) + icomp)=
     &             whole_vector_sph(idim_mid_S_stress+1,m)
               else if(icomp.eq.3.and.ir.ge.ir_CMB) then
                 c_input( ncomp * (i-1) + icomp)=
     &             whole_vector_tor(idim_mid_T_stress,m)
               end if
!the most down grid point
            else
               r_input(i)=grid_r(ir_down_stress)
               if(icomp.eq.1) then
                 c_input( ncomp * (i-1) + icomp)=
     &             whole_vector_sph(idim_down_S_stress,m)
               else if(icomp.eq.2) then
                 c_input( ncomp * (i-1) + icomp)=
     &             whole_vector_sph(idim_down_S_stress+1,m)
               else if(icomp.eq.3.and.ir.ge.ir_CMB) then
                 c_input( ncomp * (i-1) + icomp)=
     &             whole_vector_tor(idim_down_T_stress,m)
               end if
            end if
          end do
         end do
         call interpolate(ncomp,0,r_station,
     &                 r_input,c_input,c_out_interpo)
         call interpolate(ncomp,1,r_station,
     &                 r_input,c_input,dcdr_out_interpo)



cWENBO
         c1=c_out_interpo(1)
         c2=c_out_interpo(2)
         c3=c_out_interpo(3)


c**************************************************
c computing displacement
         ur=c1*Ylm(m)
         utheta=(c2*dYdtheta(m)+c3*dYdphi(m)/y)/dble(lsq)
         uphi=  (-c3*dYdtheta(m)+c2*dYdphi(m)/y)/dble(lsq)

c computing strain
         dc1dr=dcdr_out_interpo(1)
         dc2dr=dcdr_out_interpo(2)
         dc3dr=dcdr_out_interpo(3)
         du1dr=dc1dr*Ylm(m)
         du1dtheta=c1*dYdtheta(m)
         du1dphi=c1*dYdphi(m)
         du2dr=(dc2dr*dYdtheta(m)+dc3dr*dYdphi(m)/y)/dble(lsq)
         du2dtheta=(c2*d2Yd2theta(m)+c3*dphi(m)*
     &                  (-x/(y*y)*Ylm(m)+dYdtheta(m)/y))/dble(lsq)
         du2dphi=utheta*dphi(m)
         du3dr= ( -dc3dr*dYdtheta(m)+dc2dr*dYdphi(m)/y)/dble(lsq)
         du3dtheta=( -c3*d2Yd2theta(m) +
     &         c2*dphi(m)*(-x/(y*y)*Ylm(m)+dYdtheta(m)/y ) )/dble(lsq)
         du3dphi=dphi(m)*uphi

         
         epsilon11_lm=du1dr
         epsilon22_lm=(du2dtheta+ur)/grid_r(ir)

         epsilon33_lm=(du3dphi/y+utheta*x/y+ur)/grid_r(ir)
         epsilon12_lm=0.5d0*(du2dr+(du1dtheta-utheta)/grid_r(ir))
         epsilon13_lm=0.5d0*(du3dr+(du1dphi/y-uphi)/grid_r(ir))
         epsilon23_lm=0.5d0*((du3dtheta-uphi*x/y+du2dphi/y)/grid_r(ir))
!***********************Top or bottome boundary***********
!calculate the coefficients for stresses of the top and/or bottom boundary
        end if

        coef_c_lm(1)=c1
        coef_c_lm(2)=c2
        coef_c_lm(3)=c3
        coef_dcdr_lm(1)=dc1dr
        coef_dcdr_lm(2)=dc2dr
        coef_dcdr_lm(3)=dc3dr

!for convergence test
!WENBO print out for testing (r/r_source)**l decay
!        if(idep.eq.1.and.(l.eq.12000.or.l.eq.4000.or.l.eq.8000
!     &     .or.l.eq.16000.or.l.eq.1000).and.m.eq.0.and.ifreq.eq.50) then
!          do i=ngrid_r-500,ngrid_r,2
!          do i=54227,54527,3
!           print *,'coefficient=',
!     & abs(whole_vector_sph(idim_ir_sph(i),0)),l,grid_r(i),
!     &    idim_ir_sph(i)
!          end do
!        end if

        end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine comp_disp_strain0(omega,l,m,theta,phi,
     &           ir,ngrid_r,grid_r,ifreq,ir_up_stress,
     &           ir_mid_stress,ir_down_stress,
     &           ir_up_pressure,ir_mid_pressure,ir_down_pressure,
     &           ir_top_stress,ir_bot_stress,
     &           ir_top_pressure,ir_bot_pressure,
     &           ir_CMB,ir_ICB,whole_vector_sph,
     &           maxngrid_r,fluid,idim_ir_sph0,idep,idep_new,ndep_solid,
     &           ndep_fluid,grid_rho,r_station,
     &           ur,utheta,uphi,pressure_lm,epsilon11_lm,epsilon22_lm,
     &           epsilon33_lm,epsilon12_lm,epsilon13_lm,epsilon23_lm,
     &           coef_c_lm,coef_dcdr_lm,top_fluid,bot_fluid)

c variables for input/output
        implicit none
        integer ifreq
        complex*16 omega
        integer idep,idep_new,ndep_solid,ndep_fluid
        integer l,m,ir,ngrid_r,maxngrid_r,ir_CMB,ir_ICB
        integer ir_up_stress,ir_mid_stress,ir_down_stress
        integer top_fluid,bot_fluid
        integer ir_up_pressure,ir_mid_pressure,ir_down_pressure
        integer  ir_top_stress(3),ir_bot_stress(3)
        integer  ir_top_pressure(3),ir_bot_pressure(3)

        logical fluid
        real*8 theta,phi,r_station
        real*8 grid_r(*)
        real*8 grid_rho(2,maxngrid_r)
        complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
        integer idim_ir_sph0(maxngrid_r)
        
        complex*16 ur,utheta,uphi
        complex*16 pressure_lm
        complex*16 epsilon11_lm,epsilon22_lm,epsilon33_lm
        complex*16 epsilon12_lm,epsilon13_lm,epsilon23_lm
        complex*16 coef_c_lm(3),coef_dcdr_lm(3)

c other variables
        real*8 r_input(3)
        complex*16 c_input(9),c_out_interpo(3)
        complex*16 dcdr_out_interpo(3)
        real*8 factor,x,y
        complex*16 c1,dc1dr
        complex*16 dc1
        real*8 dr
        real*8 plm,plgndr
        complex*16 Ylm
        complex*16 du1dr
        integer idim_sph
        integer idim_sph_mid,idim1
        integer idim_up_pressure,idim_down_pressure
        integer idim_mid_pressure

        integer idim_up_S_stress,idim_down_S_stress
        integer idim_mid_S_stress
        integer i,icomp
        complex*16 dQdr,dQdtheta,dQdphi
        real*8 rho
        real*8 pi
        parameter ( pi=3.1415926535897932d0 )
        integer ncomp
        parameter ( ncomp= 3)



         

        if ( l.ne.0)  call error_handling(56)
c **********************************************************************
c computing the normalization factor (including the sign)
c **********************************************************************
c
        factor = dsqrt( dble(2*l+1)/(4.d0*pi) )
c **********************************************************************
c computing the displacement and its derivation
c **********************************************************************
        x=dcos(theta)
        y=dsin(theta)
        plm=plgndr(l,m,x)
        Ylm=factor*plm
c***********************************************************************
c computing the displacemetn
c***********************************************************************

c***********************************************************************
c computing the strain
c***********************************************************************
        if(fluid) then

         if(ir.eq.ir_CMB) then
             idim_sph=idim_ir_sph0(ir-1)
             rho=grid_rho(2,ir-1)
         else if(ir.eq.ir_CMB-1) then
             idim_sph=idim_ir_sph0(ir)-1
             rho=grid_rho(1,ir)
         else if(ir.lt.ir_CMB-1.and.ir.gt.ir_ICB) then
!Central differencing scheme
             idim_sph=idim_ir_sph0(ir)+1
             rho=grid_rho(1,ir)
         else if(ir.eq.ir_ICB) then
             idim_sph=idim_ir_sph0(ir)
             rho=grid_rho(1,ir)
         else
             stop 'Error in finding idim_sph in the outer core'
         end if

!idim for ir_up_pressure
         if(ir_up_pressure.eq.ir_CMB) then
             idim_up_pressure=idim_ir_sph0(ir_up_pressure-1)
         else if(ir_up_pressure.eq.ir_CMB-1) then
             idim_up_pressure=idim_ir_sph0(ir_up_pressure)-1
         else if(ir_up_pressure.lt.ir_CMB-1.and.
     &           ir_up_pressure.ge.ir_ICB) then
             idim_up_pressure=idim_ir_sph0(ir_up_pressure)
         else
            stop 'Error in finding idim_sph_up0 in the outer core'
         end if


         if(ir_mid_pressure.eq.ir_CMB) then
             idim_mid_pressure=idim_ir_sph0(ir_mid_pressure-1)
         else if(ir_mid_pressure.eq.ir_CMB-1) then
             idim_mid_pressure=idim_ir_sph0(ir_mid_pressure)-1
         else if(ir_mid_pressure.lt.ir_CMB-1.and.
     &           ir_mid_pressure.ge.ir_ICB) then
             idim_mid_pressure=idim_ir_sph0(ir_mid_pressure)
         else
            stop 'Error in finding idim_sph_mid0 in the outer core'
         end if

         if(ir_down_pressure.eq.ir_CMB) then
             idim_down_pressure=idim_ir_sph0(ir_down_pressure-1)
         else if(ir_down_pressure.eq.ir_CMB-1) then
             idim_down_pressure=idim_ir_sph0(ir_down_pressure)-1
         else if(ir_down_pressure.lt.ir_CMB-1.and.
     &           ir_down_pressure.ge.ir_ICB) then
             idim_down_pressure=idim_ir_sph0(ir_down_pressure)
         else
            stop 'Error in finding idim_sph_down0 in the outer core'
         end if

!         if(ir_up_pressure.eq.ir_CMB) then
!           rho_up=grid_rho(2,ir_up_pressure-1)
!         else
!           rho_up=grid_rho(1,ir_up_pressure)
!         end if

!         if(ir_mid_pressure.eq.ir_CMB) then
!           rho_mid=grid_rho(2,ir_mid_pressure-1)
!         else
!           rho_mid=grid_rho(1,ir_mid_pressure)
!         end if


!         if(ir_down_pressure.eq.ir_CMB) then
!           rho_down=grid_rho(2,ir_down_pressure-1)
!         else
!           rho_down=grid_rho(1,ir_down_pressure)
!         end if


         c1=whole_vector_sph(idim_sph,m)
         pressure_lm=omega*c1*Ylm
         dc1=whole_vector_sph(idim_up_pressure,m)-
     &      whole_vector_sph(idim_down_pressure,m)
         dr=grid_r(ir_up_pressure)-grid_r(ir_down_pressure)
         dc1dr=dc1/dr
         dQdr=dc1*Ylm/dr
         dQdtheta=dcmplx(0.d0,0.d0)
         dQdphi=dcmplx(0.d0,0.d0)
         ur=-dQdr/(rho*omega)
         utheta=-dQdtheta/(rho*omega*grid_r(ir))
c         uphi=-dQdphi/(rho*omega*grid_r(ir)*sin(theta))
         uphi=dcmplx(0.d0,0.d0)
c          uphi=omega*c1*Ylm



c  ***********************solid media***********************************
        else 
!for quadra inperpolation
         if(ir.eq.ir_ICB) then
             idim_sph=idim_ir_sph0(ir-1)
         else if(ir.eq.ir_ICB-1) then
             idim_sph=idim_ir_sph0(ir)-1
ccccc    idim_ir_sph(ir-1) is equal to idim_ir_sph(ir)-4
         else if(ir.gt.3) then
             idim_sph=idim_ir_sph0(ir)
         else
c             print *,'check_ir_idep',ir,idep_new
             stop 'Error in finding targeted idim_ir'
         end if


         if(ir_up_stress.eq.ir_ICB) then
            idim_up_S_stress=idim_ir_sph0(ir_up_stress-1)
         else if(ir_up_stress.eq.ir_ICB-1) then
            idim_up_S_stress=idim_ir_sph0(ir_up_stress)-1
         else if(ir_up_stress.gt.3) then
            idim_up_S_stress=idim_ir_sph0(ir_up_stress)
         else
            stop 'Error in finding targeted idim_up_stress(l=0)'
         end if

         if(ir_mid_stress.eq.ir_ICB) then
            idim_mid_S_stress=idim_ir_sph0(ir_mid_stress-1)
         else if(ir_mid_stress.eq.ir_ICB-1) then
            idim_mid_S_stress=idim_ir_sph0(ir_mid_stress)-1
         else if(ir_mid_stress.gt.3) then
            idim_mid_S_stress=idim_ir_sph0(ir_mid_stress)
         else
            stop 'Error in finding targeted idim_mid_stress'
         end if


         if(ir_down_stress.eq.ir_ICB) then
            idim_down_S_stress=idim_ir_sph0(ir_down_stress-1)
         else if(ir_down_stress.eq.ir_ICB-1) then
            idim_down_S_stress=idim_ir_sph0(ir_down_stress)-1
         else if(ir_down_stress.gt.3) then
            idim_down_S_stress=idim_ir_sph0(ir_down_stress)
         else
c            print *,'ir_down_stress',ir_down_stress,idep_new
            stop 'Error in finding targeted idim_down_stress(l=0)'
         end if

         c_input=dcmplx(0.d0,0.d0)
         do i=1,3
!The most up grid point between ir,ir_up,ir_down
          do icomp=1,1 ! only c1 is non-zero
            if(i.eq.1) then
                 r_input(i)=grid_r(ir_up_stress)
!c1
                 c_input( ncomp * (i-1) + icomp)=
     &               whole_vector_sph(idim_up_S_stress,m)
!this middle grid point
            else if(i.eq.2) then
               r_input(i)=grid_r(ir_mid_stress)
               c_input( ncomp * (i-1) + icomp)=
     &               whole_vector_sph(idim_mid_S_stress,m)
!the most down grid point
            else
               r_input(i)=grid_r(ir_down_stress)
               c_input( ncomp * (i-1) + icomp)=
     &             whole_vector_sph(idim_down_S_stress,m)
            end if
          end do
         end do
         call interpolate(ncomp,0,r_station,
     &                 r_input,c_input,c_out_interpo)
         call interpolate(ncomp,1,r_station,
     &                 r_input,c_input,dcdr_out_interpo)



cWENBO
         c1=c_out_interpo(1)


c***********************************************************************
c computing displacement and strain
c***********************************************************************
cWENBO

        ur=c1*Ylm
        utheta=dcmplx(0.d0,0.d0)
        uphi= dcmplx(0.d0,0.d0)

        dc1dr=dcdr_out_interpo(1)
        du1dr=dc1dr*Ylm
        
        epsilon11_lm=du1dr
        epsilon22_lm=ur/grid_r(ir)
        epsilon33_lm=ur/grid_r(ir)
        epsilon12_lm=dcmplx(0.d0,0.d0)
        epsilon13_lm=dcmplx(0.d0,0.d0)
        epsilon23_lm=dcmplx(0.d0,0.d0)
        end if

!***********************Top or bottome boundary***********
!calculate the coefficients for stresses of the top and/or bottom
!boundary
        coef_c_lm(1)=c1
        coef_c_lm(2)=dcmplx(0.d0,0.d0)
        coef_c_lm(3)=dcmplx(0.d0,0.d0)
        coef_dcdr_lm(1)=dc1dr
        coef_dcdr_lm(2)=dcmplx(0.d0,0.d0)
        coef_dcdr_lm(3)=dcmplx(0.d0,0.d0)

        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





        subroutine sum_disp_strain(station_disp_solid,
     &           station_disp_fluid,fluid,ur,utheta,
     &           uphi,istation_solid,istation_fluid,epsilon11,
     &           epsilon22,epsilon33,epsilon12,epsilon13,epsilon23,
     &           epsilon11_lm,epsilon22_lm,epsilon33_lm,epsilon12_lm,
     &           epsilon13_lm,epsilon23_lm,pressure,pressure_lm,
     &           copy_top_coef,copy_bot_coef,coef_c_lm,coef_dcdr_lm,
     &           l,m,nl_check_amp,i_significance,lmax,coef_c,coef_dcdr,
     &           nl_for_average)

        implicit none
        complex*16 station_disp_fluid(*)
        complex*16 station_disp_solid(*)
        integer istation_solid,istation_fluid
        logical fluid
        complex*16 ur,utheta,uphi
        complex*16 epsilon11(*),epsilon22(*),epsilon33(*)
        complex*16 epsilon12(*),epsilon13(*),epsilon23(*)
        complex*16 epsilon11_lm,epsilon22_lm,epsilon33_lm
        complex*16 epsilon12_lm,epsilon13_lm,epsilon23_lm
        complex*16 pressure(*)
        complex*16 pressure_lm
        integer lmax,l,m
        integer nl_check_amp,i_significance
        complex*16 coef_c(0:lmax,-2:2,2,3),coef_dcdr(0:lmax,-2:2,2,3)
        complex*16 coef_c_lm(3),coef_dcdr_lm(3)
        logical copy_top_coef,copy_bot_coef
        integer nl_for_average
        integer l_incheck_amp

c other variables
        real*8  fac_foraverage


        if(i_significance.eq.1) then
           fac_foraverage=1.0
        else if(i_significance.eq.0) then
          l_incheck_amp=mod(l,nl_check_amp)
          if(l_incheck_amp.eq.0) l_incheck_amp=nl_check_amp
       
          if(l_incheck_amp.le.nl_check_amp-nl_check_amp) then
             fac_foraverage=1.0
          else
             fac_foraverage=(nl_check_amp-l_incheck_amp+1.0) /
     &          nl_for_average              
          end if
        else
           write(6,*) 'Error of i_significance 
     &             (it should be either 0 or 1 here)!!'
           stop
        end if

          if(fluid) then
            station_disp_fluid(1)=
     &        station_disp_fluid(1)+ur*fac_foraverage
            station_disp_fluid(2)=
     &         station_disp_fluid(2)+utheta*fac_foraverage
           station_disp_fluid(3)=
     &         station_disp_fluid(3)+uphi*fac_foraverage

           pressure(istation_fluid)=pressure(istation_fluid) +
     &         pressure_lm*fac_foraverage

          else
            station_disp_solid(1)=
     &        station_disp_solid(1)+ur*fac_foraverage
            station_disp_solid(2)=
     &         station_disp_solid(2)+utheta*fac_foraverage
            station_disp_solid(3)=
     &         station_disp_solid(3)+uphi*fac_foraverage
            epsilon11(istation_solid)=epsilon11(istation_solid)+
     &                        epsilon11_lm*fac_foraverage
            epsilon22(istation_solid)=epsilon22(istation_solid)+
     &                        epsilon22_lm*fac_foraverage
            epsilon33(istation_solid)=epsilon33(istation_solid)+
     &                        epsilon33_lm*fac_foraverage
            epsilon12(istation_solid)=epsilon12(istation_solid)+
     &                        epsilon12_lm*fac_foraverage
            epsilon13(istation_solid)=epsilon13(istation_solid)+
     &                        epsilon13_lm*fac_foraverage
            epsilon23(istation_solid)=epsilon23(istation_solid)+
     &                        epsilon23_lm*fac_foraverage
          end if
          if(copy_top_coef) then
            coef_c(l,m,1,:)=coef_c_lm(:)
            coef_dcdr(l,m,1,:)=coef_dcdr_lm(:)
          end if
          if(copy_bot_coef) then
            coef_c(l,m,2,:)=coef_c_lm(:)
            coef_dcdr(l,m,2,:)=coef_dcdr_lm(:)
          end if



        end

	subroutine comp_vecsph_station
     &	  ( maxlmax,lmax,n_station,station_theta,station_phi,
     &	    vecsph_sph1,vecsph_sph2,vecsph_tor )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
	integer maxlmax,lmax,n_station
	real*8 station_theta(*),station_phi(*)
	complex*16 vecsph_sph1(3,0:maxlmax,-2:2,*)
	complex*16 vecsph_sph2(3,0:maxlmax,-2:2,*)
	complex*16 vecsph_tor(3,0:maxlmax,-2:2,*)
c other variables
	integer i_station
c
	if ( lmax.gt.maxlmax ) call error_handling(40)
c
	do 100 i_station=1,n_station
	  call  comp_vecsph_all
     &	  ( maxlmax,lmax,
     &	    station_theta(i_station),station_phi(i_station),
     &	    vecsph_sph1(1,0,-2,i_station),
     &	    vecsph_sph2(1,0,-2,i_station),
     &	    vecsph_tor(1,0,-2,i_station) )
  100	continue
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine comp_vecsph_all
     &	  (maxlmax,lmax,theta,phi,vecsph_sph1,vecsph_sph2,vecsph_tor)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing the vector spherical harmonics.
c   required subroutines: error_handling
c   required functions: plgndr
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
	integer maxlmax,lmax
	real*8 theta,phi
	complex*16 vecsph_sph1(3,0:maxlmax,-2:2)
	complex*16 vecsph_sph2(3,0:maxlmax,-2:2)
	complex*16 vecsph_tor(3,0:maxlmax,-2:2)
c other variables
	integer i,l,m,m0
	real*8 factor,x,plgndr
	complex*16 expimp
c constant
	real*8 pi
	parameter ( pi=3.1415926535897932d0 )
c
	do 300 l=0,lmax
	do 200 m=max0(-l,-2),min0(l,2)
c **********************************************************************
c checking the argumants
c **********************************************************************
	m0 = iabs(m)
	if ( (l.lt.0).or.(m0.gt.l) ) call error_handling(41)
c **********************************************************************
c computing the normalization factor (including the sign)
c **********************************************************************
c
	factor = 1.d0
	do 100 i=l-m0+1,l+m0
	  factor = factor * dble(i)
  100	continue
	factor = dsqrt( dble(2*l+1)/(4.d0*pi) / factor )
	if ( ( m0.ne.m ).and.( mod(m0,2).eq.1 ) ) factor = - factor
c **********************************************************************
c computing each component of the vector spherical harmonics
c **********************************************************************
	x = dcos(theta)
	expimp = cdexp( dcmplx( 0.d0, dble(m)*phi ) )
c
	vecsph_sph1(1,l,m) = factor
     &	                     * plgndr(l,m0,dcos(theta))
     &	                     * expimp
	vecsph_sph1(2,l,m) = 0.d0
	vecsph_sph1(3,l,m) = 0.d0
c
	vecsph_sph2(1,l,m) = 0.d0
	if ( l.ge.m0+1 ) then
	  vecsph_sph2(2,l,m) = factor
     &	                       * (   dble(m0) * x / dsin(theta)
     &	                             * plgndr(l,m0,x)
     &	                          + plgndr(l,m0+1,x) )
     &	                       * expimp
	else
	  vecsph_sph2(2,l,m) = factor
     &	                       * (   dble(m0) * x / dsin(theta)
     &	                             * plgndr(l,m0,x) )
     &	                       * expimp
	endif
	vecsph_sph2(3,l,m) = factor
     &	                     * dcmplx( 0.d0, dble(m) ) / dsin(theta)
     &	                     * plgndr(l,m0,dcos(theta))
     &	                     * expimp
c
	vecsph_tor(1,l,m) = dcmplx(0.d0)
	vecsph_tor(2,l,m) = vecsph_sph2(3,l,m)
	vecsph_tor(3,l,m) = - vecsph_sph2(2,l,m)
c
  200	continue
  300	continue
c	vecsph_tor(2) = factor
c     &	                * dcmplx( 0.d0, dble(m0) ) / dsin(theta)
c     &	                * plgndr(l,m0,dcos(theta))
c     &	                * expimp
c	if ( l.ge.m0+1 ) then
c	  vecsph_tor(3) = - factor
c     &	                  * (   dble(m0) * x / dsin(theta)
c     &	                        * plgndr(l,m0,x)
c     &	                      + plgndr(l,m0+1,x) )
c     &	                  * expimp
c	else
c	  vecsph_tor(3) = - factor
c     &	                  * (   dble(m) * x / dsin(theta)
c     &	                        * plgndr(l,m0,x) )
c     &	                  * expimp
c	endif
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8 function plgndr(l,m,x)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing the associated Legendre polynominal
c   (from Numerical Recipies).
c   required subroutines: error_handling
c   required functions: plgndr
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer l,m
	real*8 x
	integer i,ll
	real*8 fact,pll,pmm,pmmp1,somx2
c
	if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.) call error_handling(42)
	pmm=1.d0
	if(m.gt.0) then
	  somx2=dsqrt((1.d0-x)*(1.d0+x))
	  fact=1.d0
	  do 11 i=1,m
	    pmm=-pmm*fact*somx2
	    fact=fact+2.d0
11	  continue
	endif
	if(l.eq.m) then
	  plgndr=pmm
	else
	  pmmp1=x*dble(2*m+1)*pmm
	  if(l.eq.m+1) then
	    plgndr=pmmp1
	  else
	    do 12 ll=m+2,l
	      pll=(x*dble(2*ll-1)*pmmp1-dble(ll+m-1)*pmm)/dble(ll-m)
	      pmm=pmmp1
	      pmmp1=pll
12	    continue
	    plgndr=pll
	  endif
	endif
c
	return
	end
