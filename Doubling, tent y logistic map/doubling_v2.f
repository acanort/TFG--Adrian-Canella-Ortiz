      program doubling_v2
      implicit none

      integer i,j,t,cont_alfa
      real*8 alfa_max,alfa_min
      
      integer n_clones,t_max,num_alfa
      parameter(n_clones=20000, t_max=20000, num_alfa = 50)
      
      integer poblacion(0:t_max) !Guarda cuanta población hay en cada paso

      real*8 dt,alfa,eps,k,pi
      
      real*8 eta(n_clones)
      
      real*8 q(n_clones),p(n_clones)
      real*8 q_comp(n_clones),p_comp(n_clones)
      real*8 delta_q(n_clones),delta_p(n_clones)
      real*8 q_aux(n_clones),p_aux(n_clones)
      real*8 sep_0,sep(n_clones),sep_aux(n_clones)      

      integer*8 tau,lista_copias(10*n_clones)

!     Hacia arriba están definidas las variables que hacen correr el programa.
!     De aquí en adelante están las variables para cosas secundarias como el
!     cálculo de la presión topológica, el exponente de Lyapunov, etc.

      real*8 exp_lyap(n_clones),R(t_max),Z,exp_aux(n_clones),mu
      real*8 exp_avg, obs(n_clones), obs_avg, obs_aux(n_clones)
      
      integer i_dran,aux_int,aux_cont(0:t_max),cont
      real*8 dran_g,dran_u,aux,aux1
      call dran_ini(95537)

      open(11,file="alfa_mu_expavg_obsavg.dat")
      
      pi = 4.d0*dATAN(1.d0)   

      !Parámetros a variar: intensidad del ruido, paso temporal y exponente
      eps  = 1.d-10
      dt   = 1.d0
      k    = 1.d0
      sep_0= 1.d-3

      alfa_min =-2.5d0
      alfa_max = 2.5d0

      do cont_alfa = 0,num_alfa
         alfa = alfa_min+dfloat(cont_alfa)*(alfa_max-alfa_min)
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
         
         cont = 0
         aux_int = 0
         aux_cont = 0
         aux = 0.d0
         aux1 = 0.d0

         obs = 0.d0
         R = 0.d0
         Z = 1.d0
         mu = 0.d0
         exp_avg = 0.d0
         obs_avg = 0.d0
         
         poblacion(0) = n_clones
         
         do i=1,n_clones
            q(i) = dran_u()            
            delta_q(i) = sep_0
            q_comp(i) = q(i) + delta_q(i)                            
         enddo

         do t=1,t_max
!     write(*,*) t
            do i=1,n_clones
               eta(i) = dran_g()
               
               q(i) = MOD(2.d0*q(i),1.d0)+dsqrt(eps)*eta(i)
               
               q_comp(i) = MOD(2.d0*q_comp(i),1.d0)+dsqrt(eps)*eta(i)
               
               delta_q(i) = q_comp(i)-q(i)
               
               sep(i)=delta_q(i)/sep_0

               obs(i) = obs(i) + q(i)
               
               if (dabs(sep(i)).gt.1.d-30) then
                  exp_lyap(i)=(dlog(dabs(sep(i)))+dfloat(t-1)
     &                 *exp_lyap(i))/dfloat(t)            
               endif
               
               delta_q(i) = delta_q(i)/sep(i)
            enddo
            
            do i=1,poblacion(t-1)
               lista_copias(i) = 0
            enddo
            
            cont = 1
            do i=1,n_clones
               aux = dran_u()
               tau = floor(aux+dexp(q(i))**alfa)
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
               write(*,*) 'La poblacion ha llegado a 0.'
               stop
            endif
            
            do i=1,n_clones
               q_aux(i) = q(i)
               exp_aux(i) = exp_lyap(i)
               obs_aux(i) = obs(i)
            enddo
            
            if (poblacion(t).lt.n_clones) then 
               do i=1,poblacion(t)
                  q(i) = q_aux(lista_copias(i))
                  q_comp(i) = q_aux(lista_copias(i)) + delta_q(i)
                  
                  exp_lyap(i) = exp_aux(lista_copias(i))
                  obs(i) = obs_aux(lista_copias(i))
               enddo               
               do i=poblacion(t)+1,n_clones
                  aux_int = i_dran(poblacion(t))
                  
                  q(i) = q_aux(lista_copias(aux_int))
                  q_comp(i) = q_aux(lista_copias(aux_int)) + delta_q(i)

                  exp_lyap(i) = exp_aux(lista_copias(aux_int))
                  obs(i) = obs_aux(lista_copias(aux_int))
               enddo
            else if (poblacion(t).eq.n_clones) then
               do i=1,poblacion(t)
                  q(i) = q_aux(lista_copias(i))
                  q_comp(i) = q_aux(lista_copias(i)) + delta_q(i)

                  exp_lyap(i) = exp_aux(lista_copias(i))
                  obs(i) = obs_aux(lista_copias(i))
               enddo
            else                !if (poblacion.gt.n_clones) then
               do i=1,n_clones
                  aux_int = i_dran(poblacion(t))
                  
                  q(i) = q_aux(lista_copias(aux_int))
                  q_comp(i) = q_aux(lista_copias(aux_int)) + delta_q(i)

                  exp_lyap(i) = exp_aux(lista_copias(aux_int))
                  obs(i) = obs_aux(lista_copias(aux_int))
               enddo   
            endif

            do i=1,n_clones
               if (q(i).ge.0) then
                  q(i) = MOD(q(i),1.d0)
               else
                  q(i) = MOD(q(i),1.d0)+1.d0
               endif
            enddo
            
            R(t) = dfloat(poblacion(t))/dfloat(n_clones)
            
            Z = z*R(t)
            mu = mu+dlog(R(t))
!            write(*,*) "R",R(t),"mu",mu
         enddo                  !Fin de t

         
         
         mu = mu/dfloat(t)
         do i=1,n_clones
            obs_avg = obs_avg + obs(i)
            exp_avg = exp_avg + exp_lyap(i)
         enddo
         
         obs_avg = obs_avg/dfloat(n_clones*t_max)
         exp_avg = exp_avg/dfloat(n_clones)

         write(*,*)  alfa,mu,exp_avg,obs_avg
         write(11,*) alfa,mu,exp_avg,obs_avg                     
      enddo
      close(11)
      
      stop
      end

      
      
      include 'dranxor2.f'
