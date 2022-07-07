      program DNPendulumRK
      implicit none

      integer i,j,iter,n
      real*8 pi
      real*8 x,y,h,t,theta0,v_ang0,theta,v_ang,a_ang,l,g,alfa
      real*8 q,F_d,Omega_d,f(2),f0(2),k(4,2)

      i = 0                     !Contador para el tiempo
      h = 0.04d0                !Paso temporal
      t = 0.d0                  !Tiempo inicial
      theta = 0.d0
      g = 9.8d0
      l = 9.8d0
      n = 0

      k = 0.d0
      f = 0.d0
      f0= 0.d0
      
      pi = 4.d0*dATAN(1.d0)
      
      !Introducimos aquí las condiciones iniciales
      
      theta0 = 0.2d0             !En radianes
      v_ang0 = 0.d0
      iter = 1000000
      q = 0.5d0
      F_D = 1.2d0
      Omega_D = 2.d0/3.d0
      
      open(10,file='datospendulo.dat')
      open(11,file='seccionpoincare.dat')
      
      write(10,*) '              Tiempo',
     &            '                    Angulo',
     &            '                Posicion X',
     &            '                Posicion Y',
     &            '                 Vel. ang.'
      write(*,*)  '                      '

      
      x = l*dSIN(theta0)
      y = -l*dCOS(theta0)
      

      write(10,*) t,theta0,x,y,v_ang0,0.d0,0.d0

      f0(1) = v_ang0
      f0(2) = theta0
      
      do i=1,iter
         
         x = l*dSIN(f0(2))
         y = -l*dCOS(f0(2))

         k(1,1)= h*(-g/l*dSIN(f0(2))-q*f0(1)+F_D*dSIN(Omega_D*t))
         k(1,2)= h*(f0(1))

         k(2,1)= h*(-g/l*dSIN(f0(2)+k(1,2)/2.d0) - q*(f0(1)+k(1,1)/2.d0)
     &           + F_D*dSIN(Omega_D*(t+h/2.d0)))
         k(2,2)= h*(f0(1)+k(1,1)/2.d0)

         k(3,1)= h*(-g/l*dSIN(f0(2)+k(2,2)/2.d0) - q*(f0(1)+k(2,1)/2.d0)
     &           + F_D*dSIN(Omega_D*(t+h/2.d0))) 
         k(3,2)= h*(f0(1)+k(2,1)/2.d0)

         k(4,1)= h*(-g/l*dSIN(f0(2)+k(3,2)) - q*(f0(1)+k(3,1))
     &           + F_D*dSIN(Omega_D*(t+h))) 
         k(4,2)= h*(f0(1)+k(3,1))

         f(1) = f0(1)+1.d0/6.d0*(k(1,1)+2.d0*k(2,1)+2.d0*k(3,1)+k(4,1)) 
         f(2) = f0(2)+1.d0/6.d0*(k(1,2)+2.d0*k(2,2)+2.d0*k(3,2)+k(4,2)) 

         t = t+h
         
         if(f(2).ge.pi) then
            f(2) = f(2)-2.d0*pi
         else if(f(2).lt.-pi) then
            f(2) = f(2)+2.d0*pi
         endif

         if (t.le.60.d0) then
            write(10,*) t,f(2),x,y,f(1),0.d0,0.d0
         endif
         
         if(dabs(Omega_D*t-dfloat(n)*pi).le.(Omega_D*h)) then
            write(11,*) t,f(2),f(1)
            write(*,*) n
            n=n+1
         endif
         
         f0(1) = f(1)
         f0(2) = f(2)
         
      enddo
      
      close(10)
      close(11)
      
      stop
      end
