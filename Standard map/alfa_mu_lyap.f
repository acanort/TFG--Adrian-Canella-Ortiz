      program clones_def
      implicit none

      integer i,j,t,cont_alfa
      real*8 alfa_max,alfa_min
      
      integer n_clones,t_max,num_alfa
      parameter(n_clones=5000, t_max=500000, num_alfa=100)
      
      integer poblacion(0:t_max) !Guarda cuanta población hay en cada paso

      real*8 dt,alfa,eps,k,pi
      
      real*8 eta(n_clones)
      
      real*8 q(n_clones),p(n_clones)
      real*8 q_comp(n_clones),p_comp(n_clones)
      real*8 delta_q(n_clones),delta_p(n_clones)
      real*8 q_aux(n_clones),p_aux(n_clones)
      real*8 d_q_aux(n_clones),d_p_aux(n_clones)
      real*8 sep_0,sep(n_clones),sep_aux(n_clones)      

      integer*8 tau,lista_copias(10*n_clones)

!     Hacia arriba están definidas las variables que hacen correr el programa.
!     De aquí en adelante están las variables para cosas secundarias como el
!     cálculo de la presión topológica, el exponente de Lyapunov, etc.

      real*8 R(t_max),Z(0:num_alfa),mu(0:num_alfa)      
      real*8 exp_lyap(n_clones),exp_aux(n_clones)
      real*8 exp_avg(0:num_alfa)
      
      integer i_dran,aux_int,aux_cont(0:t_max),cont
      real*8 dran_g,dran_u,aux,aux1
      call dran_ini(123456789)

      open(11,file="alfa-mu-exp_1_1.dat")
      
      pi = 4.d0*dATAN(1.d0)      
      mu = 0.d0
      Z = 1.d0
      exp_avg = 0.d0

      !Parámetros a variar: intensidad del ruido, paso temporal y exponente
      eps  = 1.d-16
      dt   = 1.d0
      k    = 1.d0
      sep_0= 1.d-6
      alfa_max = 1.d0
      alfa_min =-1.d0

      do cont_alfa=0,num_alfa
        alfa = alfa_min + dfloat(cont_alfa)*(alfa_max-alfa_min)
     &        /dfloat(num_alfa)
        
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
         R = 0.d0
         exp_lyap = 0.d0         
         
         cont = 0
         aux_int = 0
         aux_cont = 0
         aux = 0.d0
         aux1 = 0.d0
         
         poblacion(0) = n_clones

         do i=1,n_clones
            q(i) = dran_u()
            p(i) = dran_u()-0.5d0      
            aux  = dran_u()            
            delta_q(i) = sep_0*dcos(2.d0*pi*aux) 
            delta_p(i) = sep_0*dsin(2.d0*pi*aux) 
            q_comp(i) = q(i) + delta_q(i)
            p_comp(i) = p(i) + delta_p(i)            

            exp_lyap(i) = 0.d0
         enddo

         do t=1,t_max
            do i=1,n_clones
               eta(i) = dran_g()
          
               p(i) = p(i)-dt*k*dsin(2.d0*pi*q(i))/(2.d0*pi)+dSQRT(eps)
     &              *eta(i)
               q(i) = q(i)+dt*p(i)
                        
               p_comp(i) = p_comp(i)-dt*k*dsin(2.d0*pi*q_comp(i))/(2.d0
     &              *pi) + dSQRT(eps)*eta(i)
               q_comp(i) = q_comp(i)+dt*p_comp(i)
            
               delta_q(i) = q_comp(i)-q(i)
               delta_p(i) = p_comp(i)-p(i) 

               sep(i)=dSQRT(delta_q(i)**2.d0+delta_p(i)**2.d0)/sep_0
               
               delta_q(i) = delta_q(i)/sep(i)
               delta_p(i) = delta_p(i)/sep(i)

!               if ((t.ge.0.2d0*t_max).and.(t.lt.0.8d0*t_max)) then
                  exp_lyap(i)=(dlog(abs(sep(i)))+dfloat(t-1)*
     &                 exp_lyap(i))/dfloat(t)
!               endif
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

!            write(*,*) t,poblacion(t)

            if (poblacion(t).eq.0) then
               write(*,*) 'La poblacion ha llegado a 0 en el paso',t,'.'
               stop
            endif

            do i=1,n_clones
               q_aux(i) = q(i)
               p_aux(i) = p(i)
               d_q_aux(i) = delta_q(i)
               d_p_aux(i) = delta_p(i)
               exp_aux(i) = exp_lyap(i)
            enddo

            if (poblacion(t).lt.n_clones) then 
            do i=1,poblacion(t)
               q(i) = q_aux(lista_copias(i))
               p(i) = p_aux(lista_copias(i))
               q_comp(i)=q_aux(lista_copias(i))+d_q_aux(lista_copias(i))
               p_comp(i)=p_aux(lista_copias(i))+d_p_aux(lista_copias(i))
               exp_lyap(i) = exp_aux(lista_copias(i))
            enddo               
            do i=poblacion(t)+1,n_clones
               aux_int = i_dran(poblacion(t))
               
               q(i) = q_aux(lista_copias(aux_int))
               p(i) = p_aux(lista_copias(aux_int))
               q_comp(i) = q_aux(lista_copias(aux_int))
     &              + d_q_aux(lista_copias(aux_int))
               p_comp(i) = p_aux(lista_copias(aux_int))
     &              + d_p_aux(lista_copias(aux_int))
               exp_lyap(i) = exp_aux(lista_copias(aux_int))
            enddo
         else if (poblacion(t).eq.n_clones) then
            do i=1,poblacion(t)
               q(i) = q_aux(lista_copias(i))
               p(i) = p_aux(lista_copias(i))
               q_comp(i)=q_aux(lista_copias(i))+d_q_aux(lista_copias(i))
               p_comp(i)=p_aux(lista_copias(i))+d_p_aux(lista_copias(i))
               exp_lyap(i) = exp_aux(lista_copias(i))
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
               exp_lyap(i) = exp_aux(lista_copias(aux_int))
            enddo   
         endif 

            R(t) = dfloat(poblacion(t))/dfloat(n_clones)
            Z(cont_alfa) = Z(cont_alfa)*R(t)
            mu(cont_alfa) = mu(cont_alfa) + dlog(R(t))
            

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

!               q(i) = q_aux(i) !Si lo pongo, el programa deja de funcionar :')
!               p(i) = p_aux(i)
            enddo            
         enddo                  !Fin de t

         do i=1,n_clones
            exp_avg(cont_alfa) = exp_avg(cont_alfa) + exp_lyap(i)
         enddo   
         exp_avg(cont_alfa) = exp_avg(cont_alfa)/dfloat(n_clones)
         mu(cont_alfa) = mu(cont_alfa)/(dfloat(t_max)*dt)

c$$$         aux = 0.d0
c$$$         if (cont_alfa.gt.0) then !Derivada numérica de mu
c$$$            aux = (mu(cont_alfa)-mu(cont_alfa-1))/((alfa_max-alfa_min)
c$$$     &           /dfloat(num_alfa))
c$$$         endif

         write(*,*)  alfa,mu(cont_alfa),exp_avg(cont_alfa)
         write(11,*) alfa,mu(cont_alfa),exp_avg(cont_alfa)
      enddo !Fin de cont_alfa

      close(11)
      
      stop
      end

      
      
      include 'dranxor2_new.f'
