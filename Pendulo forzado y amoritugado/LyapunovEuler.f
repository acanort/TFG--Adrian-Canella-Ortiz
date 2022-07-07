      program LyapunovEuler    !Revisar que aun no funciona :D
      implicit none

      integer i,j,iter,n,exp
      real*8 pi
      real*8 x(2),y(2),h,t,l,g,theta_ini
      real*8 theta(2),theta0(2),v_ang(2),v_ang0(2),dif_ang(0:5000)
      real*8 q,F_d,Omega_d

      i = 0                     !Contador para el tiempo
      h = 0.04d0                !Paso temporal
      t = 0.d0                  !Tiempo inicial
      theta = 0.d0
      g = 9.8d0
      l = 9.8d0

      pi = 4.d0*dATAN(1.d0)
      
      !Introducimos aquí las condiciones iniciales

      theta_ini = 0.2d0
      
      q = 0.5d0 
      F_D = 0.5d0
      Omega_D = 2.d0/3.d0

      iter = 1500
      exp = 1
c      exp  = floor(2.d0*pi/0.001d0)/10.d0 !Número de experimentos a realizar
c      write(*,*) exp
      
      open(10,file='datospendulo.dat')
      
      write(10,*) '              Tiempo',
     &      '        Diferencia angulos' 
      write(*,*) '                      '

      do n=1,exp
         theta0(1) = theta_ini         !En radianes
         theta0(2) = theta0(1)+0.001d0
         v_ang0 = 0.d0
         
         dif_ang(0) = theta0(2)-theta0(1)
         
         do i=1,iter

            do j=1,2
               v_ang(j) = v_ang0(j)-g/l*dSIN(theta0(j))*h-q*v_ang0(j)*h+ 
     &                 F_D*dSIN(Omega_D*t)*h
               theta(j) = theta0(j)+v_ang(j)*h
            enddo

            t = t+h         

            if(n.eq.1) then
               dif_ang(i) = theta(2)-theta(1)
            else
               dif_ang(i) = (dfloat(n-1)*dif_ang(i)+theta(2)-theta(1))
     &                      /dfloat(n)
            endif
               
            do j=1,2
               theta0(j) = theta(j)
               v_ang0(j) = v_ang(j)
            enddo
         enddo

         if(MOD(n,10000).eq.0) then
            write(*,*) dfloat(n)/dfloat(exp)*100
         endif

         theta_ini = theta_ini + 0.001d0
      enddo   

      t=0.d0
      do i=1,iter
         write(10,*) t,dabs(dif_ang(i))
         t=t+h
      enddo   
  
c      write(10,*) t,dABS(dif_ang),x(1),y(1),x(2),y(2)
c     &           ,0.d0,0.d0
      
      
      
      close(10)
      
      stop
      end
