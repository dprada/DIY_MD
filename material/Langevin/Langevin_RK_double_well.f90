!!!!===========================================================
!!!!    Cadena HP con analisis de componentes principales 
!!!!
!!!!     
!!!!===========================================================

PROGRAM potencia_armonico
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!     Definición variables y      !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!        lectura de datos         !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE

INTEGER::dim,n_pasos,num_divis,num_per,num_part,np,aux,aux2,random_dim
DOUBLE PRECISION, DIMENSION(:,:,:),ALLOCATABLE::gq,gp
!(num_part,coordenada,indice)
DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE::q,p,dq_dt,dp_dt,pp,qq
!(num_part,coordenada)
DOUBLE PRECISION::eps,t,Ecin,Epot
DOUBLE PRECISION, PARAMETER::delta_t=0.01d0
DOUBLE PRECISION::x_med,x2_med,x_acumul,x2_acumul,sigma_2

!# Indices para iteraciones

INTEGER::h,h2,i,j,m,m_inc,acumul


!# Parametros para Greenside-Helfand RK 3o3s2g
 
DOUBLE PRECISION, PARAMETER::A2=0.644468d0*delta_t,          &
                           & A3=0.194450d0*delta_t,          &
                           & A4=0.161082d0*delta_t,          &
                           & B21=0.516719d0*delta_t,         &
                           & B31=-0.397300d0*delta_t,        &
                           & B32=0.427690d0*delta_t,         &
                           & B41=-1.587731d0*delta_t,        &
                           & B42=1.417263d0*delta_t,         &
                           & B43=1.170469d0*delta_t,         &
                           & C01=1.0d0,                      &
                           & C12=0.271608d0,                 &
                           & C21=0.516719d0,                 &
                           & C22=0.499720d0,                 &
                           & C31=0.030390d0,                 &
                           & C32=-0.171658d0,                & 
                           & C41=1.0d0

!# Para las variables aleatorias:

REAL::beta,KT
REAL,DIMENSION(:),ALLOCATABLE::Z
INTEGER::iseed(4)

!###______Abriendo ficheros__________________________:

OPEN (11,FILE="datos.inp",STATUS="OLD",ACTION="READ")
OPEN (13,FILE="tray.oup",STATUS="NEW",ACTION="WRITE")
!OPEN (14,FILE="sigmas.oup",STATUS="NEW",ACTION="WRITE")


!###______Asignando valor____________________________:

!# Parametros de integración

READ(11,*) n_pasos
READ(11,*) m      !contador para guardar datos
READ(11,*) m_inc
                    
!# Para la ecuación a integrar:

READ(11,*) dim        
READ(11,*) num_part
READ(11,*) eps 
ALLOCATE (gq(num_part,dim,4),gp(num_part,dim,4))
ALLOCATE (q(num_part,dim),p(num_part,dim),qq(num_part,dim),pp(num_part,dim))
ALLOCATE (Z(2*dim*num_part))

gq=0.0d0
gp=0.0d0
q=0.0d0
p=0.0d0
qq=0.0d0
pp=0.0d0

Z=0.0d0


t=0.0d0 


!# Parte estocástica

READ(11,*) KT   
beta=sqrt(2*eps*KT*delta_t)

READ(11,*) iseed(1)
iseed(2)=iseed(1)+20
iseed(3)=iseed(2)+200
iseed(4)=iseed(3)+2000

aux=dim*num_part
random_dim=2*aux



q=0.0d0
p=0.0d0

q(1,1)=0.0d0

Ecin=0.0d0
Epot=0.0d0

acumul=0

x_acumul=0.0d0
x2_acumul=0.0d0
x_med=0.0d0
x2_med=0.0d0
sigma_2=0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  3o3s2g  Runge Kutta de Greenside-Helfand   !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Para termalizar::

!!$DO i=1,1000 !RK
!!$
!!$   CALL slarnv(3,iseed,random_dim,Z)
!!$   Z=beta*Z
!!$   
!!$   aux2=aux
!!$   
!!$   DO np=1,num_part
!!$      DO j=1,dim 
!!$         
!!$         aux2=aux2+1
!!$         pp(np,j)=p(np,j)+C12*Z(aux2)     !1º paso
!!$         
!!$      END DO
!!$   END DO
!!$   
!!$   CALL equation (q,pp,eps,K,gq,gp,dim,1)        
!!$   
!!$   aux2=0
!!$   
!!$   DO np=1,num_part
!!$      DO j=1,dim
!!$         
!!$         aux2=aux2+1
!!$         qq(np,j)=q(np,j)+B21*gq(np,j,1)                 !2º paso    
!!$         pp(np,j)=p(np,j)+B21*gp(np,j,1)+C21*Z(aux2)&
!!$              &+C22*Z(aux+aux2)
!!$         
!!$      END DO
!!$   END DO
!!$   
!!$   CALL equation (qq,pp,eps,K,gq,gp,dim,2)
!!$   
!!$   aux2=0
!!$      
!!$   DO np=1,num_part
!!$      DO j=1,dim
!!$         
!!$         aux2=aux2+1
!!$         qq(np,j)=q(np,j)+B31*gq(np,j,1)+B32*gq(np,j,2)  !3º paso
!!$         pp(np,j)=p(np,j)+B31*gp(np,j,1)+B32*gp(np,j,2)+C31*Z(aux2)&
!!$              &+C32*Z(aux+aux2)
!!$         
!!$      END DO
!!$   END DO
!!$   
!!$   CALL equation (qq,pp,eps,K,gq,gp,dim,3)
!!$   
!!$   aux2=0
!!$   
!!$   DO np=1,num_part
!!$      DO j=1,dim
!!$         
!!$         aux2=aux2+1
!!$         qq(np,j)=q(np,j)+B41*gq(np,j,1)+B42*gq(np,j,2)+B43*gq(np,j,3) !4º paso
!!$         pp(np,j)=p(np,j)+B41*gp(np,j,1)+B42*gp(np,j,2)+B43*gp(np,j,3)&
!!$              &+C41*Z(aux2)
!!$         
!!$      END DO
!!$   END DO
!!$   
!!$   CALL equation (qq,pp,eps,K,gq,gp,dim,4)  !fin   
!!$   
!!$   aux2=0
!!$   
!!$   DO np=1,num_part
!!$      DO j=1,dim
!!$         
!!$         aux2=aux2+1
!!$         q(np,j)=q(np,j)+A2*gq(np,j,2)+A3*gq(np,j,3)+A4*gq(np,j,4)
!!$         p(np,j)=p(np,j)+A2*gp(np,j,2)+A3*gp(np,j,3)+A4*gp(np,j,4)&
!!$              &+C01*Z(aux2)
!!$         
!!$      END DO
!!$   END DO
!!$   
!!$   t=t+delta_t 
!!$   
!!$   !# Fin del RK pero no del DO
!!$   
!!$   
!!$   
!!$END DO !RK

!!!!!

t=0.0d0


DO i=1,n_pasos !RK

   CALL slarnv(3,iseed,random_dim,Z)
   Z=beta*Z
   
   aux2=aux
   
   DO np=1,num_part
      DO j=1,dim 
         
            aux2=aux2+1
            pp(np,j)=p(np,j)+C12*Z(aux2)     !1º paso
            
         END DO
      END DO
      
      CALL equation (q,pp,eps,gq,gp,dim,1)        
      
      aux2=0
      
      DO np=1,num_part
         DO j=1,dim
            
            aux2=aux2+1
            qq(np,j)=q(np,j)+B21*gq(np,j,1)                 !2º paso    
            pp(np,j)=p(np,j)+B21*gp(np,j,1)+C21*Z(aux2)&
                 &+C22*Z(aux+aux2)
            
         END DO
      END DO
      
      CALL equation (qq,pp,eps,gq,gp,dim,2)
      
      aux2=0
      
      DO np=1,num_part
         DO j=1,dim
            
            aux2=aux2+1
            qq(np,j)=q(np,j)+B31*gq(np,j,1)+B32*gq(np,j,2)  !3º paso
            pp(np,j)=p(np,j)+B31*gp(np,j,1)+B32*gp(np,j,2)+C31*Z(aux2)&
                 &+C32*Z(aux+aux2)
            
         END DO
      END DO
      
      CALL equation (qq,pp,eps,gq,gp,dim,3)
      
      aux2=0
      
      DO np=1,num_part
         DO j=1,dim
            
            aux2=aux2+1
            qq(np,j)=q(np,j)+B41*gq(np,j,1)+B42*gq(np,j,2)+B43*gq(np,j,3) !4º paso
            pp(np,j)=p(np,j)+B41*gp(np,j,1)+B42*gp(np,j,2)+B43*gp(np,j,3)&
                 &+C41*Z(aux2)
            
         END DO
      END DO

      CALL equation (qq,pp,eps,gq,gp,dim,4)  !fin   
      
      aux2=0
      
      DO np=1,num_part
         DO j=1,dim
            
            aux2=aux2+1
            q(np,j)=q(np,j)+A2*gq(np,j,2)+A3*gq(np,j,3)+A4*gq(np,j,4)
            p(np,j)=p(np,j)+A2*gp(np,j,2)+A3*gp(np,j,3)+A4*gp(np,j,4)&
                 &+C01*Z(aux2)
            
         END DO
      END DO
      
      t=t+delta_t 
      
      !# Fin del RK pero no del DO
      

         !acumul=acumul+1
         !x_acumul=x_acumul+q(1,1)
         !x2_acumul=x2_acumul+(q(1,1))**2
         !x_med=x_acumul/(acumul*1.0d0)
         !x2_med=x2_acumul/(acumul*1.0d0)
         !sigma_2=(x2_med)-x_med**2

      
      
      IF (i==m) THEN



         WRITE(13,*) t,q(1,1)
         !WRITE(14,*) t,sigma_2
         
         m=m+m_inc          !contador para esta sección
         
      END IF



      

   END DO !RK


!###______Recogida y salida de datos despues de RK y vuelta a c.i.______:


!Para el último paso:

!Epot=K*(q(1,1))**2
!Ecin=0.5*p(1,1)**2

!print*,'----------------------------'
!print*,'q:',q(1,1),'p:',p(1,1)



CLOSE (11,STATUS="KEEP")
!CLOSE (14,STATUS="KEEP")
CLOSE (13,STATUS="KEEP")

END PROGRAM potencia_armonico

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!   Definir funcion a integrar    !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE equation (q,p,eps,dq_dt,dp_dt,dim,l)

IMPLICIT NONE

INTEGER::i,j,np
INTEGER,INTENT(IN)::dim,l
DOUBLE PRECISION,INTENT(IN)::eps
DOUBLE PRECISION,DIMENSION(1,dim),INTENT(IN)::q,p
DOUBLE PRECISION,DIMENSION(1,dim,4),INTENT(OUT)::dq_dt,dp_dt


!###___ dq_dt _______________:

dq_dt(:,:,l)=p

!###___ dp_dt _______________:

!dp_dt(:,:,l)=-eps*p(1,1)-2.0d0*K*q(1,1)
dp_dt(:,:,l)=-eps*p(1,1)-4.0d0*q(1,1)**3+4.0d0*2.0d0*q(1,1)-1.0d0


END SUBROUTINE equation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
