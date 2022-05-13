# -*- coding: utf-8 -*-

from BibliotecaFC import FluxoCarga
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 

arquivo = open('D:\\Desktop\\Estabilidade\\EstabilidadeESS\\IEEE14.txt','r')
DadosBarra,DadosLinhas,barra,ramo,Sbase,nb,nl,Ybus,B,G,teta_graus,V,teta,Pk = FluxoCarga.FluxoCarga(arquivo)

##########################################################################################################################
#                                              INICIO DO PROGRAMA
##########################################################################################################################

Contingencia = 1 # Diz qual a barra que você quer fazer a contingencia 
linha = 1 # Define a linha que se deseja retirar no periodo pos contingencia
de = 1
para = 2 

# Variaveis para integração numerica

tfout0 =  0           # Tempo inicial da falta
tfoult =  1      # Tempo de extição da falta (em segundo)
tAfter =  4.0        # Tempo maximo de simulação pós falta
maxstep = 0.001       # Passos de integração

#--------------------------------------------------------------------------------------------------------------------------
#                                   ETAPA MATRIZ REDUZIDA PRÉ FALTA
#--------------------------------------------------------------------------------------------------------------------------


# Montagem da Matriz Yl

ng = 0  
AdimiTrans = []

for i in range(len(DadosBarra['ReatanciaTransitoria'])): # Contagem do numero de geradores
    if DadosBarra['flagG'][i] != 0:
        AdimiTrans.append(-1j*DadosBarra['ReatanciaTransitoria'][i])
        ng +=1   


Y = np.zeros((ng,ng),dtype=complex) # Matriz diagona das adimitancias transitorias dos geradores (Matriz A)
AdimiTrans = np.array(AdimiTrans)
np.fill_diagonal(Y,AdimiTrans)

#print(Y)
#print('=-'*60)
#print(type(AdimiTrans))
Sconjugado = np.zeros(nb,dtype=complex) # Potencia aparente conjugado
for i in range(len(DadosBarra['num'])):
    Sconjugado[i] = DadosBarra['Pcarga'][i] + 1j * DadosBarra['Qcarga'][i]

Sconjugado = np.conjugate(Sconjugado) # Potencia aparente conjugado
Yconst = np.zeros(len(DadosBarra['num']),dtype=complex) # Carga como adimitancia constante

for i in range(len(Yconst)):
    Yconst[i] = Sconjugado[i]/V[i]**2   
    #print(Zconst[i])

YbusLPre = np.copy(Ybus)

for i in range(len(DadosBarra['num'])):
    for j in range(len(DadosBarra['num'])):
        if i == j:
            #print('entrei')
            YbusLPre[i][j] = YbusLPre[i][j] + Yconst[i]

YbusLPre_IVg = np.copy(YbusLPre)
posiG = 0 # Posição do gerador

for i in range(len(DadosBarra['num'])):
        if DadosBarra['flagG'][i] != 0:
                posiG = i
                YbusLPre[posiG,posiG] = YbusLPre[posiG,posiG] + (- 1j*DadosBarra['ReatanciaTransitoria'][i])



menosY = np.zeros((len(Y),len(DadosBarra['num'])),dtype=complex) # Matriz B e se transpor fica matriz C

Bgeradores = []

for i in range(len(DadosBarra['num'])):
    if DadosBarra['flagG'][i] != 0:
        Bgeradores.append(DadosBarra['num'][i] - 1)
         
for i in range(len(Y)):
    for j in range(len(DadosBarra['num'])):
        if j == Bgeradores[i]:
            menosY[i,j] =  1j*DadosBarra['ReatanciaTransitoria'][j]

# Montando a matriz completa

Y1Pre = np.hstack((Y,menosY))
Y2Pre = np.hstack((menosY.T,YbusLPre))
YcompletaPre = np.vstack((Y1Pre,Y2Pre))

# Montando Yreduzida para pré falta

YbusLinvPre = np.linalg.inv(YbusLPre) #Matriz D

YredPre = Y - np.dot(np.dot(menosY,YbusLinvPre),menosY.T) # Matriz Reduzida

#-------------------------------------------------------------------------------------------------------------------
#                              ETAPA MATRIZ REDUZIDA NA FALTA
#-------------------------------------------------------------------------------------------------------------------

#Contingencia = 1 # Diz qual a barra que você quer fazer a contingencia (LEMBRAR QUE O PYTHON CONTA A PARTIR DE 0)

# Calculando a matriz YredFalta para uma dada contingencia

YbusLFalta = np.copy(YbusLPre)
menosYFalta = np.copy(menosY)

for i in range(len(DadosBarra['num'])):
    for j in range(len(DadosBarra['num'])):
            YbusLFalta[i][Contingencia-1] = 0
            YbusLFalta[Contingencia-1][j] = 0

for i in range(len(DadosBarra['num'])):
    for j in range(len(DadosBarra['num'])):
        if i == Contingencia-1 and j == Contingencia-1:
            YbusLFalta[i][j] = 1

flag = 0

for i in range(ng):
        for j in range(len(DadosBarra['num'])):
                if j == Contingencia-1:
                        if menosYFalta[i][j] != 0:
                                menosYFalta[i][j] = 0
                                flag = 1
                        if flag == 1:
                                break
        if flag == 1:
                break

# Montando a matriz completa para a falta

Y1Falta = np.hstack((Y,menosYFalta))
Y2Falta = np.hstack((menosYFalta.T,YbusLFalta))
YcompletaFalta = np.vstack((Y1Falta,Y2Falta))

# Montando Yreduzida para falta

YbusLinvFalta = np.linalg.inv(YbusLFalta) #Matriz D

YredFalta = Y - np.dot(np.dot(menosYFalta,YbusLinvFalta),menosYFalta.T) # Matriz Reduzida no momento de falta 

#------------------------------------------------------------------------------------------------------------------
#                              ETAPA MATRIZ REDUZIDA PÓS FALTA
#------------------------------------------------------------------------------------------------------------------

YbusLPos = np.copy(YbusLPre)
menosYPos = np.copy(menosY)

# Tirando a linha referente a contingencia (LEMBRANDO QUE O PYTHON COMEÇA A CONTAGEM EM 0)

#linha = 2 # Define a linha que se deseja retirar no periodo pos contingencia

parcelaRetirar = complex(1+1j)
for i in range(nl):
    if i == linha-1:
        parcelaRetirar = ramo['y'][i]  # Salva a parcela do ramo que necessita ser tirada da matriz Ybus
#print(parcelaRetirar)

for i in range(nb):                    # Retira da matriz YbusLPos a parcela referente a linha onde ocorreu a contingencia
    for j in range(nb):
        if i == Contingencia-1 and j == Contingencia-1:
            YbusLPos[i,j] = YbusLPos[i,j] - parcelaRetirar
            YbusLPos[i,para-1] = 0
            YbusLPos[para-1,j] = 0
            YbusLPos[para-1,para-1] = YbusLPos[para-1,para-1] - parcelaRetirar

# Montando a matriz completa para a pos falta

Y1Pos = np.hstack((Y,menosYPos))
Y2Pos = np.hstack((menosYPos.T,YbusLPos))
YcompletaPos = np.vstack((Y1Pos,Y2Pos))

# Montando Yreduzida para falta

YbusLinvPos = np.linalg.inv(YbusLPos) #Matriz D

YredPos = Y - np.dot(np.dot(menosYPos,YbusLinvPos),menosYPos.T) # Matriz Reduzida no momento de falta 

#------------------------------------------------------------------------------------------------------------------
#                              ACHANDO O VETOR DE CORRENTES TENSÕES INTERNAS DOS GERADORES
#------------------------------------------------------------------------------------------------------------------

# Passando as tensoes e angulos para forma retangular


V_Retangular = np.ones(len(V),dtype=complex) # Vetor de tensões complexas na forma retangula
I = np.zeros(nb,dtype=complex)               # Injeção de corrente nas barras do sistema
Eg = np.zeros(len(Y),dtype=complex)          # Tensão interna dos geradores
xlinhad = []                                 # Reatancia Transitoria dos geradores
Igerador = []                                # Injeção de corrente referente as barras com geradores conectados
Vbgerador = []                               # Tensão nas barras com geradores conectados

for i in range(len(teta)):
    V_Retangular[i] = V[i]*np.cos(teta[i]) + 1j*V[i]*np.sin(teta[i])  # V*exp(j*teta) é uma representação da forma polar

I = YbusLPre_IVg @ V_Retangular # @ faz a multiplicação de matrizes

for i in range(nb):
    if DadosBarra['flagG'][i] != 0:
        Igerador.append(complex(I[i]))
        Vbgerador.append(complex(V_Retangular[i]))
        xlinhad.append(1j*(1/DadosBarra['ReatanciaTransitoria'][i]))

Igerador = np.array(Igerador)
#print(abs(Igerador))
Vbgerador = np.array(Vbgerador)
xlinhad = np.array(xlinhad)

for i in range(len(Igerador)):
    Eg[i] = Vbgerador[i] + Igerador[i]*xlinhad[i]

#print(abs(Eg))

#------------------------------------------------------------------------------------------------------------------
#                              ACHANDO AS MATRIZES C E D PARA O MOMENTO EM FALTA
#------------------------------------------------------------------------------------------------------------------

CDfalta = np.zeros((len(Y),len(Y)),dtype=complex)

for i in range(len(Y)):
    for j in range(len(Y)):
            if i != j:
                CDfalta[i][j] = abs(Eg[i])*abs(Eg[j])*YredFalta[i][j]

Cfalta = np.imag(CDfalta)
Dfalta = np.real(CDfalta)

#------------------------------------------------------------------------------------------------------------------
#                              ACHANDO AS MATRIZES C E D PARA O MOMENTO PÓS FALTA
#------------------------------------------------------------------------------------------------------------------

CDpos = np.zeros((len(Y),len(Y)),dtype=complex)

for i in range(len(Y)):
    for j in range(len(Y)):
            if i != j:
                CDpos[i][j] = abs(Eg[i])*abs(Eg[j])*YredPos[i][j]

Cpos = np.imag(CDpos)
Dpos = np.real(CDpos)

#------------------------------------------------------------------------------------------------------------------
#                              COLOCANDO AS POTENCIAS MECANICAS EM UM VETOR 
#------------------------------------------------------------------------------------------------------------------

Pm = []   # Potencia Mecanica

for i in range(nb):
    if DadosBarra['flagG'][i] != 0:
            if DadosBarra['flagG'][i] == 2: # 1 - gerador real , 2 - compensadores sincronos
                    Pm.append(0)
            else:
                    Pm.append(Pk[i]/Sbase)

Pm = np.array(Pm)

#print(Pm)

#------------------------------------------------------------------------------------------------------------------
#                              ACHANDO NOVA CONTANTE DE INERCIA E COLOCANDO DUMPING EM UM VETOR
#------------------------------------------------------------------------------------------------------------------

Dumping = []
M = []

for i in range(nb):
        if DadosBarra['flagG'][i] !=0:
                Dumping.append(DadosBarra['Dumping'][i])
                M.append(DadosBarra['ConstInercia'][i])

Dumping = np.array(Dumping)
M = np.array(M)


#------------------------------------------------------------------------------------------------------------------
#                              MONTANDO A FUNÇÃO PARA AS EQUAÇÕES DIFERENCIAIS
#------------------------------------------------------------------------------------------------------------------

w = np.zeros(ng,dtype=float)
tetaG = np.angle(Eg)
#t = 1
x0 = np.hstack((w,np.angle(Eg)))   # Estado inicial pro momento Falta

def Eqdiferencial(t,x,Pm,Yred,C,D,Eg,M):

        dw = np.ones(ng,dtype=float) 
        dDelta = np.ones(ng,dtype=float)
        G = np.real(Yred)
        w = x[0:ng]
        delta = x[ng:len(x)]
        #print((w))
        #print((delta))
        for i in range(ng):
                somatorio = []
                for j in range(ng):
                        if i != j:
                                somatorio.append(C[i][j]*np.sin(delta[i]-delta[j]) + D[i][j]*np.cos(delta[i]-delta[j]))                              
                dw[i] = (Pm[i] - (abs(Eg[i])**2) * G[i][i] - sum(somatorio)) / M[i] 
                dDelta[i] = w[i]

        x = np.hstack((dw,dDelta))

        return x


#X = Eqdiferencial(t,x0,Pm,YredFalta,Cfalta,Dfalta,Eg,M)
#print((X))

#------------------------------------------------------------------------------------------------------------------
#                              RESOLVENDO AS EQUAÇÕES DIFERENCIAIS PARA O PERIODO DE FALTA
#------------------------------------------------------------------------------------------------------------------

fun = lambda t,x: Eqdiferencial(t,x,Pm,YredFalta,Cfalta,Dfalta,Eg,M)

xFALTA = solve_ivp(fun=fun,t_span=[tfout0,tfoult],y0=x0,max_step=maxstep)

p = len(xFALTA.y[0])

x1 = xFALTA.y[:,p-1] # estado inicial do momento pos-falta

Wf = xFALTA.y[0:ng,0:len(xFALTA.y[0])]          # Omega falta
Df = xFALTA.y[ng:2*ng,0:len(xFALTA.y[0])]       # Delta falta 
tf = xFALTA.t                                   # Tempo de falta

#------------------------------------------------------------------------------------------------------------------
#                               CALCULANDO AS REFERENCIAS DO COA PARA O PERIODO FALTA
#------------------------------------------------------------------------------------------------------------------

#print(Wf[:,0])
#print(Df[:,0])
#print(len(Wf[0,:]))

# COA PARA O PERIODO FALTA:

dDeltaCOAF = []
dDeltaOmegaCOAF = []
DeltaReferenciaF = np.zeros((ng,len(Df[0,:])),dtype=float)
OmegaReferenciaF = np.zeros((ng,len(Wf[0,:])),dtype=float)

MT = sum(M)
somatorio1 = 0
somatorio2 = 0

for i in range(len(Wf[0,:])):

        for j in range(ng):
                somatorio1 = somatorio1 + M[j]*Wf[j,i]
                somatorio2 = somatorio2 + M[j]*Df[j,i]

        dDeltaCOAF.append((1/MT)*somatorio2)       # Lista contendo os valores de delta0 para cada iteração
        dDeltaOmegaCOAF.append((1/MT)*somatorio1)  # Lista contendo os valores de Omega0(W0) para cada iteração
        
        for k in range(ng):
                DeltaReferenciaF[k,i] = Df[k,i] - dDeltaCOAF[i]        # novo delta
                OmegaReferenciaF[k,i] = Wf[k,i] - dDeltaOmegaCOAF[i]   # novo omega

        somatorio1 = 0
        somatorio2 = 0

#------------------------------------------------------------------------------------------------------------------
#                                  CALCULANDO AS ENERGIAS POTENCIAIS E CINETICAS
#------------------------------------------------------------------------------------------------------------------

teta = []
for i in range(ng):
        teta.append(tetaG[i] - dDeltaCOAF[0])

somatorio3 = 0

Ec = [] # Energia Cinetica
Ep = [] # Energia Potencial

for i in range(len(Wf[0,:])):
        for j in range(ng):
                somatorio3 = somatorio3 + (1/2)*M[j] * OmegaReferenciaF[j,i]**2
        Ec.append(somatorio3)
        somatorio3 = 0

# print(Ec)

somatorio4 = 0
somatorio44 = []

for i in range(len(Wf[0,:])):
        for j in range(0,ng):
                somatorio4 = somatorio4 + ((Pm[j] - (abs(Eg[j]**2)) * np.real(YredPos[j,j])) * (DeltaReferenciaF[j,i] - teta[j]))   
        somatorio44.append(somatorio4)
        somatorio4 = 0

#print(somatorio44[0])
somatorio5 = 0
somatorio55 = 0
somatorio555 = []
#print(somatorio4)

for k in range(len(Wf[0,:])):
        for i in range(0,ng-1):
                for j in range(i+1,ng):
                        somatorio55 = somatorio55 + (Cpos[i,j] * ((np.cos(DeltaReferenciaF[i,k] - DeltaReferenciaF[j,k]) - np.cos(teta[i] - teta[j]))))
                somatorio5 = somatorio5 + somatorio55
                somatorio55 = 0
        somatorio555.append(somatorio5)
        somatorio5 = 0
        

#print(somatorio555)

somatorio6 = 0
somatorio66 = 0
somatorio667 = []

for k in range(len(Wf[0,:])):
        for i in range(0,ng-1):
                for j in range(i+1,ng):
                        if k == 0:
                                break
                        else:                                
                                somatorio66 = somatorio66 + Dpos[i,j] * ((DeltaReferenciaF[i,k]+DeltaReferenciaF[j,k]-teta[i]-teta[j]) / (DeltaReferenciaF[i,k]-DeltaReferenciaF[j,k]-teta[i]+teta[j])) * ((np.sin(DeltaReferenciaF[i,k] - DeltaReferenciaF[j,k])) - (np.sin(teta[i] - teta[j]))) 
                somatorio6 = somatorio6 + somatorio66
                somatorio66 = 0
        somatorio667.append(somatorio6)
        somatorio6 = 0


somatorio667[0] = 0
#print(len(Ec))

# Calculo da Energia Potencia

for i in range (len(Ec)):
        Ep.append(-somatorio44[i] - somatorio555[i] + somatorio667[i])


#print(V)

#------------------------------------------------------------------------------------------------------------------
#                                         CALCULANDO O VALOR DE V CRITICO
#------------------------------------------------------------------------------------------------------------------

maximo = 20
maior = 0
menor = 0
for i in range(len(Ep)):
        if Ep[i] > 0:
                maximo = Ep[i]
                if i != 0:
                        if Ep[i+1] > maximo:
                                maximo = Ep[i+1]
                                maior = Ep[i+6]
                                menor = Ep[i]

                        elif Ep[i+1] < maximo:
                                break

Vcritico = maximo
# print(Vcritico)
# print(maior)
# print(menor)

#------------------------------------------------------------------------------------------------------------------
#                                    CALCULANDO A SOMA DAS ENERGIAS E COMPARANDO COM O V CRITICO
#------------------------------------------------------------------------------------------------------------------

Et = []
t = 0
Tcritico = 0

for i in range(0,len(Ep)):
        Et.append(Ec[i]+Ep[i])

for i in range(0,len(Ep)):
        t = Et[i]
        if Et[i+1] > Vcritico:
                Tcritico = tf[i]
                break

#print(Tcritico)



#------------------------------------------------------------------------------------------------------------------
#                              RESOLVENDO AS EQUAÇÕES DIFERENCIAIS PARA O PERIODO DE PÓS FALTA
#------------------------------------------------------------------------------------------------------------------

fun1 = lambda t,x: Eqdiferencial(t,x,Pm,YredPos,Cpos,Dpos,Eg,M)

xPOS = solve_ivp(fun=fun1,t_span=[Tcritico,tAfter],y0=x1,max_step=maxstep)

p1 = len(xPOS.y[0])
#------------------------------------------------------------------------------------------------------------------
#                                         PLOTANDO OS GRAFICOS
#------------------------------------------------------------------------------------------------------------------


Wp = xPOS.y[0:ng,0:len(xPOS.y[0])]              # Omega pós
Dp = xPOS.y[ng:2*ng,0:len(xPOS.y[0])]           # Delta pós
tp = xPOS.t                                     # Tempo pós falta

#------------------------------------------------------------------------------------------------------------------
#                                    CALCULANDO AS REFERENCIAS DO COA PARA O PERIODO POS FALTA
#------------------------------------------------------------------------------------------------------------------

dDeltaCOAPF = []
dDeltaOmegaCOAPF = []
DeltaReferenciaPF = np.zeros((ng,len(Dp[0,:])),dtype=float)
OmegaReferenciaPF = np.zeros((ng,len(Wp[0,:])),dtype=float)


somatorio1 = 0
somatorio2 = 0

for i in range(len(Wp[0,:])):

        for j in range(ng):
                somatorio1 = somatorio1 + M[j]*Wp[j,i]
                somatorio2 = somatorio2 + M[j]*Dp[j,i]

        dDeltaCOAPF.append((1/MT)*somatorio2)       # Lista contendo os valores de delta0 para cada iteração
        dDeltaOmegaCOAPF.append((1/MT)*somatorio1)  # Lista contendo os valores de Omega0(W0) para cada iteração
        
        for k in range(ng):
                DeltaReferenciaPF[k,i] = Dp[k,i] - dDeltaCOAPF[i]
                OmegaReferenciaPF[k,i] = Wp[k,i] - dDeltaOmegaCOAPF[i] 

        somatorio1 = 0
        somatorio2 = 0 

pltEt = []
pltTc = []


for i in range (len(Et)):
        if Et[i] <= Vcritico:
                pltEt.append(Et[i])
                pltTc.append(tf[i])              
        elif Et[i] > Vcritico:
                break

#print(len(Ec))
#print(len(Ep))
#print(len(Et))
#print(len(tf))

#------------------------------------------------------------------------------------------------------------------
#                                         PLOTANDO OS GRAFICOS
#------------------------------------------------------------------------------------------------------------------

'''for i in range(ng):
        plt.figure(1)
        plt.plot(np.hstack((DeltaReferenciaF[i],DeltaReferenciaPF[i])),np.hstack((OmegaReferenciaF[i],OmegaReferenciaPF[i])), label=('Gerador ',i+1))  
        #plt.plot(tf,Et)   
        plt.title(' Grafico da Energia Total no tempo')
        plt.xlabel(' tempo (s)')
        plt.ylabel(' Energia (pu)')
        plt.legend()
        plt.grid(True)

        plt.figure(2)
        plt.plot(tf,np.hstack((DeltaReferenciaF[i],DeltaReferenciaPF[i])),label=('Máquina ',i+1))
        plt.title(' Grafico dos angulos no tempo ')
        plt.xlabel(' tempo (s)')
        plt.ylabel('Angulo (rad)')
        plt.legend()
        plt.grid(True)

        plt.figure(i+3)
        plt.plot(DeltaReferenciaF[i],OmegaReferenciaF[i],label=('Máquina ',i+1))
        plt.title(' Grafico de fase ')
        plt.xlabel(' Velocidade Angular (rad/s)')
        plt.ylabel('Ângulo (rad)')
        plt.legend()
        plt.grid(True)

        plt.figure(i+4)
        plt.plot(np.hstack((DeltaReferenciaF[i],DeltaReferenciaPF[i])),np.hstack((OmegaReferenciaF[i],OmegaReferenciaPF[i])),label=('Máquina ',i+1))
        plt.title(' Grafico de fase ')
        plt.xlabel(' Velocidade Angular (rad/s)')
        plt.ylabel('Ângulo (rad)')
        plt.legend()
        plt.grid(True)'''


plt.figure(1)
plt.rcParams.update({'font.size': 20})
plt.plot(pltTc,pltEt,label=('Energia Total (Et)'))
plt.plot(tf,Ec,label=('Energia Cinética (Ec)'))
plt.plot(tf,Ep,label=('Energia Potencial (Ep)'))   
plt.title(' Energia vs tepo')
plt.xlabel(' tempo (s)')
plt.ylabel(' Energia (pu)')
plt.legend()
plt.grid(True)
plt.show()



