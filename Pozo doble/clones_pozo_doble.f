      program clones_def
      implicit none

      integer i,j,cont,t
      integer n_clones,t_max
      parameter(n_clones=2000, t_max=150000)
      
      integer poblacion(0:t_max),t1,t2,t3,t4,t5,t6,t7

      real*8 dt,alfa,eps,pi
       
      real*8 eta(n_clones)
      
      real*8 q(n_clones),p(n_clones)
      real*8 q_comp(n_clones),p_comp(n_clones)
      real*8 delta_q(n_clones),delta_p(n_clones)
      real*8 sep_ini(n_clones),sep(n_clones),sep_ratio(n_clones)
      real*8 q_aux(n_clones),p_aux(n_clones),sep_aux(n_clones)
      real*8 d_q_aux(n_clones),d_p_aux(n_clones)

      integer*8 lista_copias(2*n_clones)
      
      integer i_dran,aux_int,tau
      real*8 dran_g,dran_u,aux,sep_0
      call dran_ini(957237)

      open(10,file="t1.dat")
      open(11,file="t2.dat")
      open(12,file="t3.dat")
      open(13,file="t4.dat")
      open(14,file="t5.dat")
      open(15,file="t6.dat")
      open(16,file="t7.dat")
      open(19,file="datos_gif.dat")

      t1 = 71000
      t2 = 72000
      t3 = 73000
      t4 = 74000
      t5 = 75000
      t6 = 76000

      pi = 4.d0*dATAN(1.d0)
      
      !Parámetros a variarintensidad del ruido, paso temporal y exponente
      eps  = 1.d-5
      dt   = 0.04d0
      alfa = 1.d0
      sep_0= 1.d-3

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

      poblacion(0)=n_clones
      
      do i=1,n_clones
         q(i) = -1.d0
         p(i) =  0.d0
         aux  = dran_u()
         delta_q(i) =sep_0*dcos(2.d0*pi*aux) !Hago que los clones y sus compañeros estén separados a distancia
         delta_p(i) =sep_0*dsin(2.d0*pi*aux) !sep_0 del pozo, pero en varias direcciones aleatorias
         q_comp(i) = q(i) + delta_q(i)
         p_comp(i) = p(i) + delta_p(i)
      enddo
      
      
      do t=1,t_max
         do i=1,n_clones
            !Mismo ruido para cada clon y su compañero
            eta(i) = dran_g()

            !Evoluciono cada clon y posteriormente su compañero con el mismo ruido
            q(i) = q(i) + dt*p(i)
            p(i) = p(i) + dt*4.d0*(q(i)-q(i)**3.d0)+dSQRT(eps)*eta(i)            
            
            q_comp(i) = q_comp(i) + dt*p_comp(i)
            p_comp(i) = p_comp(i) + dt*4.d0*(q_comp(i)-q_comp(i)**3.d0)
     &           +dSQRT(eps)*eta(i)

            !Caculamos la distancia entre ambos
            delta_q(i) = q_comp(i)-q(i)
            delta_p(i) = p_comp(i)-p(i) 

            !Calculamos el parámetro "separation ratio" como el cociente entre la separación del paso
            !actual y el paso anterior
            sep(i)=dSQRT(delta_q(i)**2.d0+delta_p(i)**2.d0)/sep_0

            !Normalizamos la separación entre los clones para que tenga el mismo módulo que en el paso inicial
            delta_q(i) = delta_q(i)/sep(i)
            delta_p(i) = delta_p(i)/sep(i)
         enddo

         !Aquí se copian, se eliminan o se quedan como están los clones según se indica en el paper de Nature
         !Para ello hacemos uso de un contador que recorre los elementos de un vector en el que escribiremos
         !el número de cada clon tantas veces como vayamos a copiarlo para el paso siguiente:
         !(0 veces si el clon se elimina, 1 vez si se deja como está, 2 veces si se copia)         
        
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

         if (MOD(t,floor(dfloat(t_max)/100.d0)).eq.0) then
            write(*,*)  100.d0*dfloat(t)/dfloat(t_max), poblacion(t)
         endif
      
         if (poblacion(t).eq.0) then
            write(*,*) 'La poblacion ha llegado a 0 en el paso ', t, '.'
            stop
         endif
         
         !Ahora vamos a mantener la población del sistema constante. En caso de que haya más población, eliminamos
         !de forma aleatoria copias hasta quedarnos con la población inicial. En caso de que haya menos, copiamos
         !de forma aleatoria hasta llegar a la poblacion inicial. 
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
    

         if (MOD(t,floor(dfloat(t_max)/5000.d0)).eq.0) then !Lo que hay en el denominador coincide con el numero
            do i=1,n_clones                                   !de entradas del fichero de texto para el gif
               write(19,*) t,q(i),p(i)              
            enddo
            write(19,*) ''
            write(19,*) ''
         endif

c$$$         if (dabs(dfloat(t-t1)).le.100.d0) then
c$$$            do i=1,50
c$$$               write(10,*) q(i),p(i)
c$$$            enddo   
c$$$         endif
c$$$         if (dabs(dfloat(t-t2)).le.100.d0) then 
c$$$            do i=1,50
c$$$               write(11,*) q(i),p(i)
c$$$            enddo
c$$$         endif
c$$$         if (dabs(dfloat(t-t3)).le.100.d0) then 
c$$$            do i=1,50
c$$$               write(12,*) q(i),p(i)
c$$$            enddo
c$$$         endif
c$$$         if (dabs(dfloat(t-t4)).le.100.d0) then
c$$$            do i=1,50
c$$$               write(13,*) q(i),p(i)
c$$$            enddo   
c$$$         endif
c$$$         if (dabs(dfloat(t-t5)).le.100.d0) then
c$$$            do i=1,50
c$$$               write(14,*) q(i),p(i)
c$$$            enddo   
c$$$         endif
c$$$         if (dabs(dfloat(t-t6)).le.100.d0) then
c$$$            do i=1,50
c$$$               write(15,*) q(i),p(i)
c$$$            enddo   
c$$$         endif
c$$$         if (t_max-t.lt.10000) then
c$$$            do i=1,1
c$$$               write(16,*) q(1),p(1)
c$$$            enddo
c$$$         endif
                  
      enddo     

      close(10)
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      
      stop
      end

      include 'dranxor2_new.f'
