      program clones_def
      implicit none

      integer i,j,t
      
      integer n_clones,t_max
      parameter(n_clones=1000, t_max=500000)
      
      integer poblacion(0:t_max) !Guarda cuanta población hay en cada paso

      real*8 dt,alfa,eps,k,pi
       
      real*8 eta(n_clones)
      
      real*8 q(n_clones),p(n_clones)
      real*8 q_comp(n_clones),p_comp(n_clones)
      real*8 delta_q(n_clones),delta_p(n_clones)
      real*8 q_aux(n_clones),p_aux(n_clones)
      real*8 d_q_aux(n_clones),d_p_aux(n_clones)
      real*8 sep_0,sep(n_clones),sep_aux(n_clones)

      real*8 q1(n_clones),q2(n_clones),q3(n_clones)
      real*8 q4(n_clones),q5(n_clones)
      real*8 p1(n_clones),p2(n_clones),p3(n_clones)
      real*8 p4(n_clones),p5(n_clones)

      integer*8 tau,lista_copias(10*n_clones)

!     Hacia arriba están definidas las variables que hacen correr el programa.
!     De aquí en adelante están las variables para cosas secundarias como el
!     cálculo de la presión topológica, el exponente de Lyapunov, etc.
      
      integer i_dran,aux_int,aux_cont(0:t_max),cont
      real*8 dran_g,dran_u,aux,aux1
      call dran_ini(123456789)

      open(10,file="evolucion_final.dat")
      
      pi = 4.d0*dATAN(1.d0)

      !Parámetros a variar: intensidad del ruido, paso temporal y exponente
      eps  = 1.d-16
      dt   = 0.41d0
      k    = 1.d0
      sep_0= 1.d-6
      alfa = 1.d0
      write(*,*) alfa
        
      q = 0.d0
      p = 0.d0
      q_comp = 0.d0
      p_comp = 0.d0
      delta_q = 0.d0
      delta_p = 0.d0
      q_aux = 0.d0
      p_aux = 0.d0
      eta = 0.d0
      
      sep = 0.d0
      sep_aux = 0.d0
      
      poblacion = 0         
      lista_copias = 0
      
      cont = 0
      aux_int = 0
      aux_cont = 0
      aux = 0.d0
      aux1 = 0.d0
      
      poblacion(0) = n_clones
      
      do i=1,n_clones
!         q(i) = 0.01d0*dran_u()+0.495d0
!         p(i) = 0.01d0*(dran_u()-0.5d0)
         q(i) = 0.5d0
         p(i) = 0.0d0
         aux  = dran_u()            
         delta_q(i) = sep_0*dcos(2.d0*pi*aux) 
         delta_p(i) = sep_0*dsin(2.d0*pi*aux) 
         q_comp(i) = q(i) + delta_q(i)
         p_comp(i) = p(i) + delta_p(i)                            
      enddo

      do t=1,t_max
         do i=1,n_clones
            eta(i) = dran_g()
            
            p(i) = p(i)-dt*k*dsin(2.d0*pi*q(i))/(2.d0*pi)+dSQRT(eps)
     &           *eta(i)
            q(i) = q(i)+dt*p(i)
            
            p_comp(i) = p_comp(i)-dt*k*dsin(2.d0*pi*q_comp(i))/(2.d0
     &           *pi) + dSQRT(eps)*eta(i)
            q_comp(i) = q_comp(i)+dt*p_comp(i)
            
            delta_q(i) = q_comp(i)-q(i)
            delta_p(i) = p_comp(i)-p(i) 
            
            sep(i)=dSQRT(delta_q(i)**2.d0+delta_p(i)**2.d0)/sep_0
            
            delta_q(i) = delta_q(i)/sep(i)
            delta_p(i) = delta_p(i)/sep(i)
         enddo
         
         do i=1,poblacion(t-1)
            lista_copias(i) = 0
         enddo
         
         cont = 1
         do i=1,n_clones
            aux = dran_u()
            tau = floor(aux+sep(i)**alfa)
            if (tau.gt.0) then
               do j=1,tau
                  lista_copias(cont)=i
                  cont = cont+1                  
               enddo
               poblacion(t) = poblacion(t) + tau
            endif
         enddo
         
!     write(*,*) t,poblacion(t)
         
         if (poblacion(t).eq.0) then
            write(*,*) 'La poblacion ha llegado a 0 en el paso ', t, '.'
            stop
         endif
         
         do i=1,n_clones
            q_aux(i) = q(i)
            p_aux(i) = p(i)
            d_q_aux(i) = delta_q(i)
            d_p_aux(i) = delta_p(i)
         enddo
         
         if (poblacion(t).lt.n_clones) then 
            do i=1,poblacion(t)
               q(i) = q_aux(lista_copias(i))
               p(i) = p_aux(lista_copias(i))
               q_comp(i)=q_aux(lista_copias(i))+d_q_aux(lista_copias(i))
               p_comp(i)=p_aux(lista_copias(i))+d_p_aux(lista_copias(i))
            enddo               
            do i=poblacion(t)+1,n_clones
               aux_int = i_dran(poblacion(t))
               
               q(i) = q_aux(lista_copias(aux_int))
               p(i) = p_aux(lista_copias(aux_int))
               q_comp(i) = q_aux(lista_copias(aux_int))
     &              + d_q_aux(lista_copias(aux_int))
               p_comp(i) = p_aux(lista_copias(aux_int))
     &              + d_p_aux(lista_copias(aux_int))
            enddo
         else if (poblacion(t).eq.n_clones) then
            do i=1,poblacion(t)
               q(i) = q_aux(lista_copias(i))
               p(i) = p_aux(lista_copias(i))
               q_comp(i)=q_aux(lista_copias(i))+d_q_aux(lista_copias(i))
               p_comp(i)=p_aux(lista_copias(i))+d_p_aux(lista_copias(i))
            enddo
         else                   !if (poblacion.gt.n_clones) then
            do i=1,n_clones
               aux_int = i_dran(poblacion(t))
               
               q(i) = q_aux(lista_copias(aux_int))
               p(i) = p_aux(lista_copias(aux_int))
               q_comp(i) = q_aux(lista_copias(aux_int))
     &              + d_q_aux(lista_copias(aux_int))
               p_comp(i) = p_aux(lista_copias(aux_int))
     &              + d_p_aux(lista_copias(aux_int))
            enddo   
         endif       
         
         
         do i=1,n_clones
            q_aux(i) = MOD(q(i),1.d0)
            if (q_aux(i).lt.0.d0) then
               q_aux(i) = q_aux(i) + 1.d0
            endif
            
            if (p(i)+0.5d0.ge.0.d0) then
               p_aux(i) = MOD(p(i)+0.5d0,1.d0)-0.5d0
            else
               p_aux(i) = MOD(p(i)+0.5d0,1.d0)+0.5d0
            endif
            
!     q(i) = q_aux(i) !Si lo pongo, el programa deja de funcionar :')
!     p(i) = p_aux(i)
         enddo

         if (t.ge.t_max-5000) then
            do i=1,20
               write(10,*) q_aux(i),p_aux(i)
            enddo
         endif
         
      enddo                     !Fin de t
      
      close(10)
           
      stop
      end

      
      
      include 'dranxor2_new.f'
