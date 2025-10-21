!    - Criar ficheiro K.dat contendo a matriz NxM (exemplo: 2x3)
!    - Criar ficheiro F.dat contendo o vetor Mx1 (exemplo: 2x1)
!    - Garantir que os valores sao reais (usar .e0) para evitar lixo numerico
!    - Se fossem inteiros, inicializacao nao seria necessaria
!    - Caso nao soubessemos o tipo, inicializar com 0 para evitar valores residuais

! Escrever o programa principal ('main')
!    - Declarar 'implicit none' para obrigar a declaracao de todas as variaveis
!    - Declarar variaveis:
!        K(2,3)  -> matriz de 2 linhas e 3 colunas
!        F(2)    -> vetor de 2 elementos
!        i, j    -> contadores de linhas e colunas

!Subrotina para ler a matriz K
!    - Entrada: nome do ficheiro (ficheiroK) [intent(in)]
!    - Saida: matriz K [intent(out)]
!    - Abrir o ficheiro K.dat usando unidade 10
!    - Para cada linha i=1 a N:
!        - Ler os elementos da linha i (j=1 a M) e armazenar em K(i,j)
!    - Fechar o ficheiro

! Subrotina para ler o vetor F
!    - Entrada: nome do ficheiro (ficheiroF) [intent(in)]
!    - Saida: vetor F [intent(out)]
!    - Abrir o ficheiro F.dat usando unidade 11
!    - Para cada elemento j=1 a M:
!        - Ler o elemento F(j) e armazenar no vetor
!    - Fechar o ficheiro

! Ler os ficheiros na 'main'
!    - Chamar subrotina lerK("K.dat", K)
!    - Imprimir a matriz K linha a linha usando um loop
!    - Chamar subrotina lerF("F.dat", F)
!    - Imprimir o vetor F elemento a elemento

! Comentarios adicionais
!    - Utilizar 'do' para percorrer linhas da matriz ou elementos do vetor
!    - Intent(in) serve para variaveis que sao apenas lidas
!    - Intent(out) serve para variaveis que serao preenchidas na subrotina
!    - Usar unidades de ficheiro fixas (10 e 11) para simplificar a leitura
!    - Garantir que os ficheiros existem antes de chamar as subrotinas

program EscreverMatriz
    implicit none
    real*8 :: K(2,3), F(2)         ! matriz K 2x3 e vetor F 2x1 (E SE quisessemos defenir o tamanho das matrizes )
    real*8 :: KT(3,2), FT(1:2)
    integer :: i,j                   ! contador de linhas

    ! Matriz K
    call lerK("K.dat",K)            ! chama a subrotina para ler a matriz K do ficheiro
    print*, "Matriz K:"             ! cabecalho para impressao da matriz
    do i = 1, 2
        print*, K(i,:)              ! imprime linha i completa da matriz K
    end do

    ! Vetor F
    call lerF("F.dat",F)            ! chama a subrotina para ler o vetor F do ficheiro
    print*, "Vetor F:"              ! cabecalho para impressao do vetor
    print*, F                       ! imprime todos os elementos do vetor F

    call transpostaF (F, FT)
    call transpostaK (K, KT)

    open(12, file="OUTPUT.dat", status="unknown", action="write")
!status="unknown" -tenta abrir se existir, ou cria se não exist
!tatus="new" - cria um ficheiro novo, e dá erro se o ficheiro já existir.
!status="replace" - cria um ficheiro novo, e se o ficheiro já existir ele é substituído.
!status="unknown" - tenta abrir se existir, ou cria se não exist

    write(12,*) "Matriz K transposta (3x2):"
    do i = 1, size(KT,1)
        write(12, *) KT(i,:)     ! escreve linha de KT
    end do

    write(12,*) "Vetor F transposto (1x2):"
    write(12,*) FT(1,:)

    close (12)

end program EscreverMatriz

!------------------------------------------------------------
! Subrotina para ler a matriz K
! K - output (a matriz lida)
! ficheiroK - input (nome do ficheiro)
!------------------------------------------------------------
subroutine lerK(ficheiroK,K)
    character(len=*), intent(in) :: ficheiroK   ! nome do ficheiro passado pela main (entrada !TEMOS QUE FAZER SEMPRE ISTO ?
    real*8, intent(out) :: K(2,3)              ! matriz de saida 2x3
    integer :: i,j                               ! indices de leitura

    open(10, file=ficheiroK, status="old", action="read")  ! abre o ficheiro K.dat para leitura

    do i = 1, 2                                   ! percorre as 2 linhas da matriz
        read(10,*) (K(i,j), j=1,3)               ! le os 3 elementos de cada linha
    end do

    close(10)                                     ! fecha o ficheiro
end subroutine lerK

!|-----------------------------------------------------------------------------------------------------------
!| intent (in) -- A variavel entra na subrotina apenas para ser lida. Nao pode ser modificada dentro da subrotina.
!| intent (out) -- A variavel vai ser escrita/modificada dentro da subrotina e usada no programa principal.
!|                   Quando entra na subrotina, o valor antigo e ignorado.
!| Vamos usar intent(in) para o nome do ficheiro pois vamos utilizar o seu nome na main
!| Vamos usar intent(out) para outras variaveis pois serao alteradas dentro da subrotina e depois usadas na main
!|_________________________________________________________________________________________________________

!------------------------------------------------------------
! Subrotina para ler o vetor F
! F - output (o vetor lido)
! ficheirof - input (nome do ficheiro)
!------------------------------------------------------------
subroutine lerF(ficheirof,F)
    character(len=*), intent(in) :: ficheirof   ! nome do ficheiro passado pela main
    real*8, intent(out) :: F(2)                 ! vetor de saida 2x1
    integer :: j                                 ! contador

    open(11, file=ficheirof, status="old", action="read")   ! abre o ficheiro F.dat para leitura

    do j = 1, 2                                   ! percorre os 2 elementos do vetor
        read(11,*) F(j)                           ! le o elemento j do vetor
    end do

    close(11)                                     ! fecha o ficheiro
end subroutine lerF


subroutine transpostaK (K,KT)
        real*8, intent(in) :: K(2,3)
    real*8, intent(out) :: KT(3,2)
    integer :: i,j

    do i = 1,2
        do j = 1,3
            KT(j,i) = K(i,j)
        end do
    end do

end subroutine

subroutine transpostaF (F,FT)
    real*8, intent(in) :: F(2)
    real*8, intent(out) :: FT(1,2)
    integer :: i

    do i = 1,2
        FT(1,i) = F(i)
    end do

end subroutine
