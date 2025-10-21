import numpy as np
import matplotlib.pyplot as plt

def resolver_viga():
    print("=== Análise de Viga com Cargas Distribuídas e Pontuais ===")
    L = input_positivo("Insira o comprimento da viga [m]: ")
    a = input_a(L)
    b = input_b(a, L)
    q = input_positivo("Insira o valor da carga distribuída q [kN/m]: ")
    nforcas = input_numero_forcas()

    # Cargas pontuais
    cargas_pontuais = []
    for i in range(1, nforcas + 1):
        P = input_positivo(f"Insira a intensidade da carga P_{i} [kN]: ")
        c = input_distancia_forca(L, i)
        cargas_pontuais.append((P, c))

    # Carga distribuída equivalente
    Q = q * (b - a)
    xQ = a + (b - a) / 2

    # Reações nos apoios
    Reacao_A = (sum(P * (L - c) for P, c in cargas_pontuais) + Q * (L - xQ)) / L
    Reacao_B = (sum(P for P, _ in cargas_pontuais) + Q) - Reacao_A

    print(f"\nReações de apoio: RA = {Reacao_A:.2f} kN, RB = {Reacao_B:.2f} kN")

    # Funções de esforço transverso e momento fletor
    def esforco_transverso(x):
        V = Reacao_A
        for P, c in cargas_pontuais:
            if x >= c:
                V -= P
        if a <= x <= b:
            V -= q * (x - a)
        elif x > b:
            V -= q * (b - a)
        return V

    def momento_fletor(x):
        M = Reacao_A * x
        for P, c in cargas_pontuais:
            if x >= c:
                M -= P * (x - c)
        if a <= x <= b:
            w = q * (x - a)
            centroide = (x - a) / 2
            M -= w * centroide
        elif x > b:
            w = q * (b - a)
            centroide = (b - a) / 2
            M -= w * (x - (a + centroide))
        return M

    x_vals = np.linspace(0, L, 500)
    V_vals = [esforco_transverso(x) for x in x_vals]
    M_vals = [momento_fletor(x) for x in x_vals]

    # Desenahr diagramas
    plt.figure(figsize=(12, 6))
    plt.subplot(2, 1, 1)
    plt.plot(x_vals, V_vals, 'b', label='Esforço Transverso V(x)',color="orange")
    plt.axhline(0, color='k', linewidth=0.8)
    plt.ylabel('V (kN)')
    plt.grid(True)
    plt.legend()
    plt.subplot(2, 1, 2)
    plt.plot(x_vals, M_vals, 'r', label='Momento Fletor M(x)', color="blue")
    plt.axhline(0, color='k', linewidth=0.8)
    plt.xlabel('x (m)')
    plt.ylabel('M (kN·m)')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

#Introduzir valores positivos
def input_positivo(mensagem):
    valor = float(input(mensagem))
    while valor <= 0:
        print("Valor inválido! Insira um valor positivo.")
        valor = float(input(mensagem))
    return valor

#Introduzir um valor valido entre 0 e L
def input_a(L):
    a = float(input("Insira o valor de a em metros (>= 0 e <= L): "))
    while a < 0 or a > L:
        print("Valor inválido! a deve ser >= 0 e <= L.")
        a = float(input("Insira o valor de a em metros (>= 0 e <= L): "))
    return a

#Introduzir o valor de b de maneira que seja maior que o valor anteriormente defenido (a) e menor que L
def input_b(a, L):
    b = float(input(f"Insira o valor de b em metros (>= {a} e <= L): "))
    while b < a or b > L:
        print(f"Valor inválido! b deve ser >= {a} e <= {L}.")
        b = float(input(f"Insira o valor de b em metros (>= {a} e <= L): "))
    return b

#Introduzir um numero inteiro positivo de forças
def input_numero_forcas():
    nforcas = input("Insira o número de forças: ")
    while not nforcas.isdigit() or int(nforcas) <= 0:
        print("Valor inválido! Deve ser um número inteiro positivo.")
        nforcas = input("Insira o número de forças: ")
    return int(nforcas)

#Introduzir o valor da distancia de cada força, com as devidas restriçoes
def input_distancia_forca(L, contador):
    distancia = float(input(f"Qual a distância da carga P_{contador} ao ponto A (<= {L}): "))
    while distancia > L or distancia < 0:
        print(f"Distância inválida! Deve ser entre 0 e {L}.")
        distancia = float(input(f"Qual a distância da carga P_{contador} ao ponto A (<= {L}): "))
    return distancia

if __name__ == "__main__":
    resolver_viga()