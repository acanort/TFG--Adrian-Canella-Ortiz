      program histo_lyap_def
      implicit none

      integer*8 i,j,t,t1,t2
      integer*8 t_max,it,n_cajas
      parameter(t_max=100000,it=50000,n_cajas=200)
      
      real*8 k,dt
      real*8 q,p,q_comp,p_comp,q_aux,p_aux
      real*8 delta_q, delta_p

      real*8 sep,sep_0
      real*8 exp_lyap
      real*8 exp_lyap_t1(it),exp_lyap_t2(it),exp_lyap_t3(it) !Guardamos en el tiempo t_i el exp.lyap de todas las ejecuciones del sistema

      real*8 media
      
      real*8 pi

      real*8 histo_min,histo_max,ancho
      integer*8 histo1(n_cajas),histo2(n_cajas),histo3(n_cajas)
      integer*8 sum1,sum2,sum3, cont_break
      
      real*8 dran_g,dran_u,aux
      call dran_ini(95537)

      open(11,file='n1-exp.dat')
      open(12,file='n2-exp.dat')
      open(13,file='n3-exp.dat')
      open(21,file='histo1.dat')
      open(22,file='histo2.dat')
      open(23,file='histo3.dat')
      open(30,file='exp_lyap.dat')
      
      pi = 4.d0*dATAN(1.d0)
      histo1 = 0
      histo2 = 0
      histo3 = 0

      dt   = 1.d0
      k    = 1.d0
      sep_0 = 1.d-10

      histo_max = 0.8d0
      histo_min = 0.0d0
      ancho = (histo_max - histo_min)/dfloat(n_cajas)

      write(*,*) dt,k
      
      t1 = 20
      t2 = 300

      write(*,*) 'Tiempos a los que se hará el histograma:'
      write(*,*) t1,t2,t_max
      
      

      media = 0.d0
      cont_break = 0
      
      do j=1,it         
         q = dran_u()
         p = dran_u()-0.5d0
         aux  = dran_u()
         delta_q = sep_0*dcos(2.d0*pi*aux)
         delta_p = sep_0*dsin(2.d0*pi*aux)
         q_comp = q + delta_q
         p_comp = p + delta_p

         exp_lyap = 0.d0

         do t=1,t_max
            p = p - k*dt*dsin(2.d0*pi*q)/(2.d0*pi)
            q = q + dt*p

            p_comp = p_comp - k*dt*dsin(2.d0*pi*q_comp)/(2.d0*pi)
            q_comp = q_comp + dt*p_comp

            delta_q = q - q_comp
            delta_p = p - p_comp
        
            sep = dsqrt(delta_q**2.d0 + delta_p**2.d0)/sep_0
            
            delta_q = delta_q/sep
            delta_p = delta_q/sep

c$$$            if (dabs(sep).gt.1.d-30) then
c$$$               exp_lyap = (dlog(sep)+dfloat(t-1)*exp_lyap)/dfloat(t)
c$$$            endif
            
            if (dabs(sep).gt.1.d-30) then
               exp_lyap = exp_lyap + dlog(dabs(sep))
            endif
            

            q_comp = q + delta_q
            p_comp = p + delta_p
            
            if (t.eq.t1) then
               exp_lyap_t1(j) = exp_lyap/dfloat(t1)
            endif
            if (t.eq.t2) then
               exp_lyap_t2(j) = exp_lyap/dfloat(t2)
            endif

            if (mod(j,10).eq.0) then
               write(30,*) t, exp_lyap/dfloat(t)
            endif
         enddo                  !Final de la evolucion temporal de cada clon         
         
         write(30,*) ''
         write(30,*) ''

         exp_lyap = exp_lyap/dfloat(t_max)         
         exp_lyap_t3(j) = exp_lyap
         

         write(11,*) t1,exp_lyap_t1(j)
         write(12,*) t2,exp_lyap_t2(j)
         write(13,*) t ,exp_lyap_t3(j)
        
         do i=1,n_cajas
            if(exp_lyap_t1(j).le.histo_min+dfloat(i)*ancho) then
               histo1(i) = histo1(i)+1
               exit
            endif
         enddo
         
         do i=1,n_cajas
            if(exp_lyap_t2(j).le.histo_min+dfloat(i)*ancho) then
               histo2(i) = histo2(i)+1
               exit
            endif
         enddo
         
         do i=1,n_cajas
            if(exp_lyap_t3(j).le.histo_min+dfloat(i)*ancho) then
               histo3(i) = histo3(i)+1
               exit
            endif
         enddo

         media = media + exp_lyap
         write(*,*) exp_lyap,media
      enddo !Final de cada una de las simulaciones independientes de cada clon

      write(*,*) 'cont break',cont_break
      media = media/dfloat(it-cont_break)

      write(*,*) "exp_lyap medio = ", media
      
      sum1 = 0 !Comprobamos que todos los histogramas se hayan llenado con todas los lanzamientos del sistema
      sum2 = 0
      sum3 = 0
      do i=1,n_cajas
         aux = histo_min+dfloat(i)*ancho
         write(21,*) aux,histo1(i)
         sum1 = sum1+histo1(i)
         write(22,*) aux,histo2(i)
         sum2 = sum2+histo2(i)
         write(23,*) aux,histo3(i)
         sum3 = sum3+histo3(i)
      enddo
      write(*,*) ''
      write(*,*) sum1,sum2,sum3
     
      close(11)
      close(12)
      close(13)
      close(21)
      close(22)
      close(23)
      close(30)
      
      stop
      end

      include 'dranxor2_new.f'
