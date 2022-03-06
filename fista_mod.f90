MODULE fista_mod
  
  USE precision, WP => DP
  IMPLICIT NONE
  
CONTAINS
  
  SUBROUTINE fista(A,B,s_crit,lambda,Xk)
    
    !==================================================================
    ! SUBROUTINE fista: Applies Fast Iterative Shrinkage-Threshold    | 
    !                   Algoritm to minimize:                         |
    !                                                                 |
    !                   J=||AX-B||² + Lambda||X||_1                   |
    !                                                                 | 
    ! IN:             A: Direct model matrix                          |
    !                 B: Observable                                   |
    !                 s_crit: stop criterion                          |
    !                 lambda: trade-off parameter                     |
    !                                                                 |
    !                                                                 |
    ! OUT:            Xk: Estimated sparse solution                   |
    !==================================================================
    
    REAL(KIND=WP), INTENT(IN), DIMENSION(:,:)::A
    REAL(KIND=WP), INTENT(IN), DIMENSION(:)::B
    REAL(KIND=WP), INTENT(OUT), DIMENSION(:)::Xk
    
    REAL(KIND=WP), ALLOCATABLE::AT_A(:,:),Xkl1(:),Yk(:)
    REAL(KIND=WP)::alpha,lambda,tau,tk,tkp1,tkl1,misfit_f,l1norm_f,&
         misfit_temp
    INTEGER::steps,i,j,k,s_crit
    LOGICAL::stop_criterion,file_status

    !=================================================================
    !Se verifica que Xk y A tengan dimensiones compatibles
    !=================================================================


    IF(SIZE(A,1)  /= SIZE(B) .OR. SIZE(A,2) /= SIZE(Xk))THEN 
       WRITE(*,*)'Las dimensiones de los arreglos no son compatibles para aplicar FISTA'
       STOP
    END IF

    !=================================================================
    !Se busca el maximo autovalor de AT_A.
    !El archivo "max_alpha" tendra el valor pre-calculado,
    !si este archivo no existe el autovalor se calcula ahora
    !=================================================================

    INQUIRE(FILE='max_alpha', EXIST=file_status)
    IF(file_status)THEN
       OPEN(UNIT=10,FILE='max_alpha',ACTION='read')
       READ(10,*)alpha
       CLOSE(10)
    ELSE
       WRITE(*,*)'El archivo "max_alpha" no existe, se calcula el maximo autovalor'
       ALLOCATE(AT_A(SIZE(A,2),SIZE(A,2)))
       AT_A=MATMUL(TRANSPOSE(A),A)
       steps=50
       CALL findeigen(AT_A,alpha,steps)
       DEALLOCATE(AT_A)
    END IF
    alpha=alpha*1.1_WP
    write(*,'(A20,f10.2)')'Max autovalor:',alpha

    !===FISTA=========================================================

    ALLOCATE(Xkl1(SIZE(A,2)),Yk(SIZE(A,2)))
      
    tau=lambda/(2.0_WP*alpha)
    Xk=0.0_WP
    Xkl1=0.0_WP
    Yk=0.0_WP
    tk=1.0_WP
    tkp1=0.0_WP
    tkl1=0.0_WP
    misfit_temp=10_WP**10
    k=0
    stop_criterion=.TRUE.

    !OPEN(UNIT=20,FILE='MISFIT')

    DO WHILE(stop_criterion)
       
       k=k+1

       Xk=st_function(Yk-&
            (1.0_WP/alpha)*MATMUL(TRANSPOSE(A),MATMUL(A,Yk)-B),tau)
       
       tkp1=(1.0_WP+SQRT(1.0_WP+4.0_WP*(tk**2)))/2.0_WP
       Yk=Xk+(tkl1/tkp1)*(Xk-Xkl1)
         
       Xkl1=Xk
       tkl1=tk
       tk=tkp1
       
       misfit_f=DOT_PRODUCT(MATMUL(A,Xk)-B,&
            MATMUL(A,Xk)-B)
       l1norm_f=SUM(ABS(Xk))
      
       !WRITE(20,*)misfit_f,l1norm_f

       IF(MOD(k,100)==0) WRITE(*,'(I5,5X,F10.6,5X,F10.6)')k,misfit_f,l1norm_f
       
       !criterios de corte!
       IF(s_crit < 0 .AND. MOD(k,100)==0 .AND. &
            ABS(misfit_temp-misfit_f)/misfit_f < 10.0_WP**(s_crit)) THEN
          stop_criterion=.FALSE. 
          WRITE(*,*)'Criterio de corte: misfit estable'
       ELSE
          misfit_temp=misfit_f
       END IF
       
       IF(s_crit > 0 .AND. k > s_crit) THEN
          stop_criterion=.FALSE.
          WRITE(*,*)'Criterio de corte: numero maximo de iteraciones alcanzado'
       END IF
    
    END DO

    CLOSE(20)

    DEALLOCATE(Xkl1,Yk)

  END SUBROUTINE fista
  
  !========================================================================

  SUBROUTINE findeigen(MATRIX,eigenval,steps)

    !esta subrutina encuentra el maximo autovalor de la MATRIX usando el
    !metodo de las potencias de Rayleigth.

    IMPLICIT NONE
    INTEGER, INTENT(IN)::steps  
    REAL(KIND=WP), INTENT(IN), DIMENSION(:,:)::MATRIX 
    REAL(KIND=WP), INTENT(OUT) :: eigenval 
    
    REAL(KIND=WP), ALLOCATABLE:: X(:) 
    INTEGER :: i, j

    ALLOCATE(X(SIZE(MATRIX,1)))

    X  = 1.0_WP 
    
    DO i = 1, steps
       X=MATMUL(MATRIX,X)
       eigenval=MAXVAL(X)
       IF(eigenval == 0) EXIT			
       X=X/eigenval
    END DO

  
    DEALLOCATE(X)
    
  END SUBROUTINE findeigen

  !========================================================================

  FUNCTION st_function(VECTOR,tau_in)
    
    !Shrinkage-Tresholding Function

    IMPLICIT NONE
    REAL(KIND=WP)::tau_in
    REAL(KIND=WP), DIMENSION(:)::VECTOR
    REAL(KIND=WP), ALLOCATABLE, DIMENSION(:)::st_function
    INTEGER i

    ALLOCATE(st_function(SIZE(VECTOR)))
    
    DO i=1,SIZE(VECTOR)
       st_function(i)=MAX(0.0_WP,ABS(VECTOR(i))-tau_in)*SIGN(1.0_WP,VECTOR(i))
    END DO

    RETURN

    DEALLOCATE(st_function)

  END FUNCTION st_function



END MODULE fista_mod
