module globals
	implicit none	
	integer,parameter :: N=2, f=100000, n1 = 500, fooo = 15, n_max=10, foooo = 16, fooooo = 17, filesum = 5000
	real (8) t, y(N), C(n1), n_mol, n_dot,  sum_alpha2, R_max, R(n1), n_mol_temp1, n_mol_temp2, r_temp1, r_temp2
	real (8) Pg_dot, a_n(2*n_max), gamma, T_K(n1), V, RT(n1), U_dot, U, pre_y(N)
	real(8), parameter :: dt=1.0d0*10.0d0**(-11.0d0) , Pi=3.141592650d0, Rg=8.3d0, T_K0=293.0d0, rho=1000.0d0
	real(8), parameter :: P0=100000.0d0, nu=1.0d0*10.0d0**(-6.0d0), D=1.44d0*10.0d0**(-9.0d0), gamma0=0.07275d0
	real(8), parameter :: H=0.0003d0, P_mugen=250000.0d0, Pv0=249990.0d0, A=1.2d0, dx=1.0d0*10.0d0**(-8), alpha = Pi
	real(8), parameter :: wc = 1500.0d0, CO2 = 44.0d0, dr=1.0d0*10.0d0**(-10.0d0)
	real(8), parameter :: br0 = 2 * gamma0 / (Pv0 - P0)* 1.0d0, Cv=29.0d0
	real(8), parameter :: D_T=2.6d0*10.0d0**(-5.0d0), Tdx=3.0d0*10.0d0**(-7.0d0), D_T2=2.6d0*10.0d0**(-5.0d0)
	real(8), parameter :: t_max=1.0d0/1000*10, tstart=0.0d0/f, t2start=30.0d0/f
end module globals

module subfunctions
	use globals
	implicit none
	contains
	! 
	function RP_func(t1,y1,n_mol1,Pg_dot1, T_K1) result (f1)
		real(8) t1, y1(N), f1(N), e1, e2, e3, e4, e5, n_mol1, Pg_dot1,T_K1(n1)

		e1 = Rg / rho
		e2 = P0 / rho
		e3 = 4 * nu
		e4 = 2 * gamma / rho
		e5 = 1.0d0 / rho / wc

		V = 4.0d0 * Pi * y1(1)**3 / 3

		f1(1)= y1(2)
		f1(2)= e1 * n_mol1 * T_K1(1) / V / y1(1) - e2 * (1 + A * sin(2 * Pi * f * t1 + alpha)) / y1(1) &
			- 3 * y1(2)**2 / 2 / y1(1) - e3 * y1(2) / y1(1)**2 - e4 / y1(1)**2 
	 	!y1(1)=radius  y1(2)=velocity  f(1)=rdot  f1(2)=vdot​

	end function RP_func

	function dif_func(y1,n_mol1,C1,R1,T_K1) result (f1)
		!real(8), parameter :: 
		real(8) t1, y1(N), f1(n1*2+1), C1(n1), C_plus_dt(n1), C_plus_dt2(n1), R1(n1), R_plus_dt(n1), J, n_dot1, n_mol1, alpha2
		real(8) dR(n1), u_R(n1), T_K1(n1)
		integer i
		
		C_plus_dt(1) = H * 3 * n_mol1 * Rg * T_K(1) / 4 / Pi / y1(1)**3
		C_plus_dt2(1) = H * 3 * n_mol1 * Rg * T_K(1) / 4 / Pi / y1(1)**3
		C_plus_dt(n1) = H * P_mugen
		C_plus_dt2(n1) = H * P_mugen
		R_plus_dt(1) = y1(1)

		do i = 2, n1
			u_R(i) = y1(1)**2 / R1(i)**2 * y1(2)
			R_plus_dt(i) = R1(i) + u_R(i) * dt
			C_plus_dt(i) = C1(i)
		end do
		
		
		do i = 2, n1-1
			C_plus_dt2(i) = C_plus_dt(i) + 2 * D / R_plus_dt(i) * (C_plus_dt(i) - C_plus_dt(i-1)) / (R_plus_dt(i) - R_plus_dt(i-1)) * dt &
						+ D * ((C_plus_dt(i+1) - C_plus_dt(i)) / (R_plus_dt(i+1) - R_plus_dt(i)) - &
						(C_plus_dt(i) - C_plus_dt(i-1)) / (R_plus_dt(i) - R_plus_dt(i-1))) / (R_plus_dt(i) - R_plus_dt(i-1)) * dt
		end do

		J = D * (C_plus_dt2(2) - C_plus_dt2(1)) /  (R_plus_dt(2) - R_plus_dt(1))

		n_dot1 = 4 * Pi * y1(1)**2 * J

		f1(1) = n_dot1
		do i = 2, n1+1
			f1(i) = C_plus_dt2(i-1)
		end do
		do i = n1+2, n1*2+1
			f1(i) = R_plus_dt(i-(n1+1))
		end do

	end function dif_func

	function temperature_func(y1,RT1,T_K1) result (f1)
		!real(8), parameter :: 
		real(8) t1, y1(N), f1(2*n1+1), RT1(n1), u_R(n1)
		real(8) T_K1(n1), T_plus_dt(n1), T_plus_dt2(n1), R_plus_dt(n1), J, U_dot1
		integer i

		T_plus_dt(1) = T_K1(1)
		T_plus_dt2(1) = T_K1(1)
		T_plus_dt(n1) = T_K0
		T_plus_dt2(n1) = T_K0
		
		R_plus_dt(1) = y1(1)

		do i = 2, n1
			u_R(i) = y1(1)**2 / RT1(i)**2 * y1(2)
			R_plus_dt(i) = RT1(i) + u_R(i) * dt
			T_plus_dt(i) = T_K1(i)
		end do

		do i = 2, n1-1
			T_plus_dt2(i) = T_plus_dt(i) + 2 * D_T / R_plus_dt(i) * (T_plus_dt(i) - T_plus_dt(i-1)) / (R_plus_dt(i) - R_plus_dt(i-1)) * dt &
						+ D_T * ((T_plus_dt(i+1) - T_plus_dt(i)) / (R_plus_dt(i+1) - R_plus_dt(i)) - &
						(T_plus_dt(i) - T_plus_dt(i-1)) / (R_plus_dt(i) - R_plus_dt(i-1))) / (R_plus_dt(i) - R_plus_dt(i-1)) * dt
		end do

		J = D_T2 * (T_plus_dt2(2) - T_plus_dt2(1)) /  (R_plus_dt(2) - R_plus_dt(1))

		U_dot1 = 4 * Pi * y1(1)**2 * J

		f1(1) = U_dot1

		do i = 2, n1+1
			f1(i) = T_plus_dt2(i-1)
		end do
		do i = n1+2, n1*2+1
			f1(i) = R_plus_dt(i-(n1+1))
		end do

	end function temperature_func

	function spherical_harmonics(y1,a_n1,t1,n_mol1,Pg_dot1,T_K1) result (f1)
		real(8) f1(2*n_max), y1(2), delta, a_n1(2*n_max), r2dot(2), t1, n_mol1, Pg_dot1, rhog, Anchil, Bnchil, An, Bn
		real(8) T_K1(n1)
		integer i
		
		do i = 1, n_max
			delta = min((nu / (2 * Pi * f))**(1.0d0/2), y1(1) / 2 / i)
			r2dot = RP_func(t1,y1,n_mol1,Pg_dot1,T_K1)
			f1(2*i - 1) = a_n1(2*i)
			rhog = 3 * n_mol * CO2 / 1000 / 4 / Pi / y1(1)**3
			Anchil = (i-1)*(r2dot(2)/y1(1) - (i+1)*(i+2)*gamma/rho/y1(1)**3 - 2*nu*y1(2)/y1(1)**3 * &
					((i+1)*(i+2) - i*(i+2)/(1+2*delta/y1(1))))
			Bnchil = 3 * y1(2) / y1(1) + 2 * nu / y1(1)**2 * (-(i-1)*(i+1)*(i+2) + i*(i+2)**2 / (1+2*delta/y1(1)))
			An =  rho/(i+1) / (rhog/i + rho/(i+1)) * Anchil - rhog/i / (rhog/i + rho/(i+1)) * (i-1) * y1(2)**2 / y1(1)
			Bn = rho/(i+1) / (rhog/i + rho/(i+1)) * Bnchil
			f1(2*i) =  -1.0d0 * Bn * a_n1(2*i) + An * a_n1(2*i-1)
		end do

	end function spherical_harmonics


end module subfunctions

module subroutines
	use globals
	use subfunctions
	implicit none
	contains
	
	
	subroutine Runge_Kutta_RP
		real(8) k1(N),k2(N),k3(N),k4(N)
		k1(1:N)=RP_func(t,y(1:N),n_mol,Pg_dot,T_K)
		k2(1:N)=RP_func(t+1/2*dt,y(1:N)+1/2*dt*k1(1:N),n_mol,Pg_dot,T_K)
		k3(1:N)=RP_func(t+1/2*dt,y(1:N)+1/2*dt*k2(1:N),n_mol,Pg_dot,T_K)
		k4(1:N)=RP_func(t+dt,y(1:N)+dt*k3(1:N),n_mol,Pg_dot,T_K)
		y(1:N) = y(1:N) + dt*( k1(1:N) + 2*k2(1:N) + 2*k3(1:N) + k4(1:N) )/6
	end subroutine

	subroutine Runge_Kutta_SH
		real(8) k1(2*n_max),k2(2*n_max),k3(2*n_max),k4(2*n_max)
		k1(1:2*n_max)=spherical_harmonics(y(1:N),a_n(1:2*n_max),t,n_mol,Pg_dot,T_K)
		k2(1:2*n_max)=spherical_harmonics(y(1:N), a_n(1:2*n_max)+1/2*dt*k1(1:2*n_max), t,n_mol,Pg_dot,T_K)
		k3(1:2*n_max)=spherical_harmonics(y(1:N), a_n(1:2*n_max)+1/2*dt*k2(1:2*n_max), t,n_mol,Pg_dot,T_K)
		k4(1:2*n_max)=spherical_harmonics(y(1:N), a_n(1:2*n_max)+dt*k3(1:2*n_max), t,n_mol,Pg_dot,T_K)
		a_n(1:2*n_max) = a_n(1:2*n_max) + dt*( k1(1:2*n_max) + 2*k2(1:2*n_max) + 2*k3(1:2*n_max) + k4(1:2*n_max) )/6
	end subroutine

	subroutine Euler_dif
		real(8) f(n1*2+1)
		integer i
		f = dif_func(y,n_mol,C,R,T_K)
		n_dot = f(1)
		do i = 1, n1
			C(i) = f(i+1)
		end do
		do i = 1, n1
			R(i) = f(i+n1+1)
		end do
		n_mol = n_mol + dt * n_dot
	end subroutine

	subroutine Euler_temperture
		real(8) f(n1*2+1)
		integer i

		f = temperature_func(y,RT,T_K)
		U_dot = f(1)
		do i = 2, n1
			T_K(i) = f(i+1)
		end do
		do i = 1, n1
			RT(i) = f(i+n1+1)
		end do

		U = U + dt * U_dot - n_mol * Rg * T_K(1) / (4.0d0 * Pi * y(1)**3 / 3) * (4.0d0 * Pi / 3) * (y(1)**3 - pre_y(1)**3)

		T_K(1) = U / Cv / n_mol
	end subroutine
end module subroutines

program RK
	use globals
	use subroutines
	use subfunctions
	implicit none
	integer(8) i, imax, j, k, filenum, orbitnum
	integer,parameter :: fo=11, foo=13, jmax=20, kmax=60
	real(8),parameter::  xrange=5.0d0, yrange=10.0d0*10.0d0**(-6.0d0)
	real(8) a2dot, SH(2*n_max), g1, g2, g3, g4, g5, x_vec, y_vec, xx, yy, v_plus, v_minus, rr
	real(8) orbit_r(500), orbit_v(500),TF(2*n1+1),Pv, PdV
	character(40) ffout,fffout, ffffout, fout

	filenum = 0
	orbitnum = 1
	open (11,file="diffusion_convection_Rayleigh.d")
	
	write(*,*)br0
	y(1) = br0
	y(2) = 0.0d0      !initial velocity
	T_K(:) = T_K0
	n_mol = 4 * Pi * br0**3 * Pv0 / 3 / Rg / T_K(1)	!initial n
	C(:) = H * P_mugen
	U = Cv * T_K(1) * n_mol
	
	

	do i = 1, n_max
		a_n(2*i - 1) = 0.1d0 * br0
		a_n(2*i) = 0.1d0 * (P0/rho)**(1.0d0/2)
	end do

	do i = 1, n1
		R(i) = br0 + dx * (i - 1)
	end do

	do i = 1, n1
		RT(i) = br0 + Tdx * (i - 1)
	end do
	

	imax=t_max/dt


	do i=1,imax
		t=i*dt

		if (mod(i, int(imax/0.5/10**3)) == 1)then
			write(fout,'("DCR2_",i4.4,".d")') i/int(imax/0.5/10**3)
			open(foo,file=fout)
			do j = 1, n1-2
				write(foo,*) R(j), C(j), sin(2 * Pi * f * t + alpha)
			end do
			! write(fout,'("DCT2_",i4.4,".d")') i/int(imax/0.5/10**3)
			! open(foo,file=fout)
			! do j = 1, n1
			! 	write(foo,*) RT(j), T_K(j), sin(2 * Pi * f * t + alpha)
			! end do
		end if
		if (mod(i,int((imax-int(tstart/dt))/filesum/10)) == 1 .and. t>tstart)then
			Pv = n_mol * Rg * T_K(1) / (4.0d0 / 3 * Pi * y(1)**3)
			PdV = Pv * (4.0d0 * Pi / 3) * (y(1)**3 - pre_y(1)**3)
			!write(fo,*) t, y(1),y(2),T_K(1),U_dot!n_dot/4/Pi/y(1)**2!, RP_func(t,y,n_mol,Pg_dot)!a_n!y(1),y(2),sin(2 * Pi * f * t + alpha)
			!write(fo,*)t,Pv,4.0d0 / 3 * Pi * y(1)**3,T_K(1), Pv*4.0d0 / 3 * Pi * y(1)**3,n_mol,y(1),U_dot,U_dot*dt,PdV,y(1)**3 - pre_y(1)**3
			write(fo,*)t,y(1),y(2),T_K(1),n_mol
			if(t>tstart)then
				write(ffout,'("nullcline2_",i4.4,".d")') filenum
				write(fffout,'("nullcline2_orbit_",i4.4,".d")')filenum
				write(ffffout,'("nullcline2_vector_",i4.4,".d")')filenum
				open(fooo,file=ffout)
				open(foooo,file=fffout)
				open(fooooo,file=ffffout)
				a2dot = -1.0d0 / y(1) * (3 * y(2) + 18 * nu / (y(1) + 2 * min((nu / (2 * Pi * f))**(1.0d0/2), y(1) / 2)))
				SH(1:2*n_max) = spherical_harmonics(y(1:N),a_n(1:2*n_max),t,n_mol,Pg_dot,T_K)
				
				g1 = 3 * Rg * T_K(1) / 4 / rho / Pi
				g2 = P0 / rho
				g3 = 4 * nu
				g4 = 2 * gamma0 / rho
				g5 = 1.0d0 / rho / wc

				do j = -jmax, jmax
					do k = 0, kmax
						xx = xrange * j / jmax
						yy = yrange * k / kmax
						x_vec = g1 * n_mol / yy**4 - g2 * (1 + A * sin(2 * Pi * f * t + alpha)) / yy &
							- 3 * xx**2 / 2 / yy - g3 * xx / yy**2 - g4 / yy**2
						y_vec = xx
						write(fooooo,*)xx, yy*10**6, x_vec*3*10.0d0**(-9.0d0), y_vec*10.0d0**(-7.0d0)*10**6
					end do
				end do
				orbit_r(orbitnum) = y(1)
				orbit_v(orbitnum) = y(2)
				do j = 10000, 40000
					rr = j * dr
					v_plus = - 4 * nu / 3 / rr + (16 * nu**2 / 9 / rr**2 + Rg * T_K(1) / 2 / rho / Pi * n_mol / rr**3 &
						- 2 * P0 / 3 / rho * (1 + A * sin(2 * Pi * f * t + alpha)) -  4 * gamma0 / 3 / rho / rr)**(1.0d0/2)
					v_minus = - 4 * nu / 3 / rr - (16 * nu**2 / 9 / rr**2 + Rg * T_K(1) / 2 / rho / Pi * n_mol / rr**3 &
						- 2 * P0 / 3 / rho * (1 + A * sin(2 * Pi * f * t + alpha)) -  4 * gamma0 / 3 / rho / rr)**(1.0d0/2)
					write(fooo,*)v_plus, v_minus, rr*10**6, y(2), y(1)**10**6, xrange*10**6/30
				end do
				filenum = filenum + 1
				orbitnum = orbitnum + 1

				do j = 1, orbitnum-1
					write(foooo,*)orbit_r(j)*10**6, orbit_v(j)
				end do
					


			end if
			
		 end if
		pre_y = y
		call Runge_Kutta_RP
		call Euler_dif
		call Runge_Kutta_SH
		call Euler_temperture
		
		
		gamma = -0.00015392d0 * T_K(1) + 0.117666d0

		if ( y(1)<0 )then
			y(1)=0.0d0
			write(*,*) 'r<0'
			write(*,*)t*f,'周期',t,'秒で終わる'
			exit
		end if
    
	end do
	write(*,*)t*f,'周期',t,'秒'
	close(fo)
	
end program RK



