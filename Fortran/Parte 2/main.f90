




!para nao usar o module,nem interface,nem contais a solucao e criar uma subrotina apenas para ler as dimensoes
program main
    implicit none
    integer :: NL,NC
    real*8, allocatable :: A(:,:),B(:,:),R(:,:),RT(:,:),X(:,:)
    !MATRIZ A
    call LerA("A.dat",A,NL,NC)
    print*, "A matriz A:"
    call EscreverA (A,NL,NC)
    !MATRIZ B
    call LerB("B.dat",B,NL,NC)
    print*,"A matriz B:"
        call EscreverB (B,NL,NC)

    call Cholesky(A,R)
    call TransCholesky(R,RT)

    call Resolve (R,RT,B,X)
    call Guardar ("resolve.dat",A,B,R,RT,X,NL,NC)

        deallocate (A,B,R,RT,X)
    contains

subroutine LerA(ficheiroA,m,l,c)
    implicit none
!começar ler o numero de linhas e colunas
    character (len=*),intent(in) ::ficheiroA
    real*8, allocatable:: m(:,:)
    integer:: l,c,i,j
        open(10,file=ficheiroA,status="old",action="read")
            read (10,*) L, C
            allocate(m(l,c))
            do i=1,L
                read(10,*)(m(i,j),j=1,C)
            end do
        Close(10)
end subroutine

subroutine EscreverA (m,l,c)
    implicit none
    real*8,intent(in):: m(l,c)
    integer :: l,c,i,j
        do i=1,l
                write(*,*)(m(i,j),j=1,c)
        end do
end subroutine

subroutine LerB (ficheiroB,v,linhas,colunas)
    implicit none
!começar ler o numero de linhas e colunas
    character (len=*),intent(in) :: ficheiroB
    real*8, allocatable :: v(:,:)
    integer :: linhas,colunas,i,j
        open(11,file=ficheiroB,status="old",action="read")
            read(11,*) linhas,colunas
            allocate(v(linhas,colunas))
            do i=1,linhas
                read(11,*)(v(i,j),j=1,colunas)
            end do
        close(11)
end subroutine

subroutine EscreverB(v,linhas,colunas)
    implicit none
    real*8 :: v(linhas,colunas)
    integer :: linhas,colunas,i,j
        do i=1,linhas
                write(*,*)(V(i,j),j=1,colunas)
        end do
end subroutine

subroutine Cholesky (A,R)
    implicit none
    real*8, intent(in):: A(:,:)
    real*8, allocatable,intent(out):: R(:,:)
    integer ::i,j,k
    real*8 :: sum
    sum=0.e0 !inicializar a soma em 0 para nao haver lixo numerico
    R=0.e0 !Inicializar a matriz R toda em 0
    !agora vamos criar uma rotina por R para conseguir escrever as matrizes

    do i=1,size(A,1)  !percorrer as linhas
        do j=1,i !percorrer as colunas
            do k=1,i-1
                sum=sum+R(k,i)*R(k,j)
            end do
            if(i==j) then !diagonal principal
            R(i,i)=sqrt(A(i,i)-sum)
            else
                R(i,j)=(A(i,i)-sum)/R(j,j)
            end if
        end do
    end do
end subroutine

subroutine TransCholesky (R,RT)
    implicit none
    real*8, intent(in)::R(:,:)
    real*8, intent(out):: RT(:,:)
    integer ::i,j
        do i = 1,size(R,1)
            do j = 1,size (R,1)
            RT(j,i) = R(i,j)
        end do
    end do
end subroutine

subroutine Resolve (R,RT,B,X)
implicit none
    !-------------------------------------------------------------
    ! Resolve o sistema A*X = B, onde A = RT * R (Cholesky)
    !
    ! Entradas:
    !   R   - matriz triangular superior da decomposição
    !   RT  - transposta de R (triangular inferior)
    !   B   - termos independentes
    !
    ! Saída:
    !   X   - solução do sistema
    !-------------------------------------------------------------

    real*8, intent(in)  :: R(:,:), RT(:,:), B(:,:)
    real*8, intent(out) :: X(:,:)
    integer :: n, m
    integer :: i, j, k
    real*8, allocatable :: Y(:,:)
    real*8 :: soma

    n = size(R,1)        ! dimensão do sistema
    m = size(B,2)        ! número de colunas em B (podem ser vários sistemas)

    allocate(Y(n,m))     ! vetor/matriz auxiliar para a solução intermédia
    X = 0.0d0
    Y = 0.0d0

    ! Etapa 1: resolver RT * Y = B (substituição direta)
    do k = 1, m
        do i = 1, n
            soma = B(i,k)
            do j = 1, i-1
                soma = soma - RT(i,j) * Y(j,k)
            end do
            Y(i,k) = soma / RT(i,i)
        end do
    end do

    ! Etapa 2: resolver R * X = Y (retro-substituição)
    do k = 1, m
        do i = n, 1, -1
            soma = Y(i,k)
            do j = i+1, n
                soma = soma - R(i,j) * X(j,k)
            end do
            X(i,k) = soma / R(i,i)
        end do
    end do

    deallocate(Y)
end subroutine

subroutine Guardar(Resolve,A,B,R,RT,X,l,c)
    implicit none
    character(len=*) :: Resolve
    integer :: l,c,i,j
    real*8 :: A(:,:), B(:,:),R(:,:),RT(:,:),X(:,:)
    open(12,file="Resolve.dat",status="replace",action="write")
    write(12,*) '--- Matriz A ---'
    do i = 1, l
        write(12,'(100(f12.6,1x))') (A(i,j), j = 1, l)
    end do

    write(12,*)
    write(12,*) '--- Matriz B ---'
    do i = 1, l
        write(12,'(100(f12.6,1x))') (B(i,j), j = 1, c)
    end do

    write(12,*)
    write(12,*) '--- Matriz R (Cholesky, triangular superior) ---'
    do i = 1, l
        write(12,'(100(f12.6,1x))') (R(i,j), j = 1, l)
    end do

    write(12,*)
    write(12,*) '--- Matriz RT (transposta de R) ---'
    do i = 1, l
        write(12,'(100(f12.6,1x))') (RT(i,j), j = 1, l)
    end do

    write(12,*)
    write(12,*) '--- Solucao X do sistema A * X = B ---'
    do i = 1, l
        write(12,'(100(f12.6,1x))') (X(i,j), j = 1, c)
    end do

    close(12)
 end subroutine

end program
