# -*- coding: utf-8 -*-

import numpy as np

class LeitorDados: 
    def ler_dados(arq):
        #arq = open('C:\\Users\\joao\\Desktop\\Mestrado\\disciplinas_mestrado\\Estabilidade em SEP\\Trabalho2\\IEEE14.txt', 'r')
        texto = arq.readlines()
        arq.close()
        b = []
        c = []
        aux = 0
        Sbase = 100
        for i in range(len(texto)):
            if texto[i] == '9999\n':
                aux += 1
                continue
            if aux == 0:                                   #------ Dados de Barra ------- append - add um elemento no final da lista
                a = []
                a.append(texto[i][0:4])                    # N da barra
                a.append(texto[i][7:8])                    # Tipo da barra (2 - ref; 1 - PV; 0-PQ)
                a.append(texto[i][22:26])                  # Tensao Inicial
                a.append(texto[i][27:30])                  # Angulo inicial
                a.append(texto[i][30:35])                  # Potencia Ativa Gerada
                a.append(texto[i][35:40])                  # Potencia Reativa Gerada
                a.append(texto[i][40:46])                  # Potencia Reativa minima
                a.append(texto[i][46:51])                  # Potencia Reativa maxima
                a.append(texto[i][56:60])                  # Potencia Ativa consumida
                a.append(texto[i][60:65])                  # Potencia Reativa consumida
                a.append(texto[i][65:70])                  # Shunt da barra
                a.append(texto[i][70:74])                  # flag indicando que há gerador na barra
                a.append(texto[i][74:80])                  # constante de inercia H do gerador
                a.append(texto[i][80:82])                  # Dumping referente ao gerador
                a.append(texto[i][82:91])                  # Reatancia transitoria do gerador
                b.append(a[:])
            if aux == 1:                                   # -------- Dados das linhas ------ 
                a = []
                a.append(texto[i][0:4])                    # De barra x 
                a.append(texto[i][4:12])                   # Para barra y
                a.append(texto[i][16:23])                  # Resistencia da linha [%]
                a.append(texto[i][23:29])                  # Reatancia da linha  [%]
                a.append(texto[i][29:35])                  # Reantancia Shunt da linha [%]
                a.append(texto[i][35:40])                  # Tap dos transformadores [%]
                a.append(texto[i][40:45])                  # Tap minimo
                a.append(texto[i][45:50])                  # Tap maximo
                a.append(texto[i][50:55])                  # Angulo trafo defasador
                a.append(texto[i][55:64])                  # Normal MVA rating
                a.append(texto[i][64:68])                  # Valor MVA emergencial
                c.append(a[:])

        for j in range(len(b)):                            # Fazer o fatiamento e eliminando espacos em branco dos dados das barras
            for k in range(len(b[j])):
                b[j][k] = b[j][k].split()

        for j in range(len(c)):                            # Fazer o fatiamento e eliminando espacos em branco dos dados das linhas
            for k in range(len(c[j])):
                c[j][k] = c[j][k].split()

        for i in range(len(b)):                            # Transformando as variaveis em tipo float e substituindo espacos vazios por 0.0
            for j in range(len(b[i])):
                for k in range(len(b[i][j])):
                    b[i][j][k] = float(b[i][j][k])
                if b[i][j] == []:
                    b[i][j] = [0.0] 
                

        for i in range(len(c)):                            # Substituo espaços vazios por 0
            for j in range(len(c[i])):
                for k in range(len(c[i][j])):          
                    c[i][j][k] = float(c[i][j][k])
                if  c[i][j] == []:                
                    c[i][j] = [0.0]            


        # Cria um dicionario organizando todas as informacoes contidas na lista b das barras:

        DadosBarra = {'num': None,'tipo':None,'tensao':None,'angulo':None,'Pger':None,'Qger': None,'Qmin_ger':None,'Qmax_ger':None,'Pcarga':None,'Qcarga':None,'Qshunt': None,'flagG': None,'ConstInercia':None,'Dumping':None,'ReatanciaTransitoria':None}
        cont1 = 0
        for k in DadosBarra.keys():
            temp = []
            for i in range(len(b)):
                temp.append(b[i][cont1][0])
            cont1 += 1
            DadosBarra[k] =  temp[:]
            del(temp)

        for i in range(len(DadosBarra['ReatanciaTransitoria'])):
            if DadosBarra['ReatanciaTransitoria'][i] != 0:
                DadosBarra['ReatanciaTransitoria'][i] = 1/DadosBarra['ReatanciaTransitoria'][i]

        #print(DadosBarra['flagG'])
        #print('-='*100)
        #print(DadosBarra['ConstInercia'])
        #print('-='*100)
        #print(DadosBarra['Dumping'])
        #print('-='*100)
        #print(DadosBarra['ReatanciaTransitoria'])
        # Cria um dicioanrio organizando todas as informacoes contidas na lista c das linhas:

        DadosLinhas = {'De': None,'Para':None,'R': None,'Xl': None,'Xc': None,'tap': None,'tapmin':None,'tapmax':None,'fi':None}
        cont2 = 0
        for k in DadosLinhas.keys():
            temp = []
            for i in range(len(c)):
                temp.append(c[i][cont2][0])
            cont2 += 1
            DadosLinhas[k] = temp[:]
            del(temp)

        for k in range(len(DadosLinhas['De'])):                        #Passando numeros das linhas para serem iguais ao enderecamento de memoria do python
            DadosLinhas['De'][k] -= 1
        for k in range(len(DadosLinhas['Para'])):
            DadosLinhas['Para'][k] -= 1

        #print(DadosLinhas['R'])
        #print('=='*80)
        #print(DadosLinhas['Xl'])
        #print('=='*80)
        #print(DadosLinhas['tap'])
        #print(DadosBarra['Qshunt'])
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        nb = len(DadosBarra['num'])  # Obtem o numero de barras do sistema
        nl = len(DadosLinhas['De'])  # Obtem o numero de linhas do sistema
        npv = 0
        npq = 0
        for i in range(nb):
            if DadosBarra['tipo'] == 1:
                npv += 1
            elif DadosBarra['tipo'] == 0:
                npq += 1 

        #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        # Ajustando os dados de linha e de barra para o uso correto no programa
        # Calculo de P e Q liquido

        temp3 = []
        temp4 = []
        barra = {'P':None,'Q':None,'bsh':None}
        for i in range(nb):
            temp3.append(((DadosBarra['Pger'][i])-(DadosBarra['Pcarga'][i]))/Sbase)
            temp4.append(((DadosBarra['Qger'][i])-(DadosBarra['Qcarga'][i]))/Sbase)
            DadosBarra['Pcarga'][i] = DadosBarra['Pcarga'][i]/Sbase
            DadosBarra['Qcarga'][i] = DadosBarra['Qcarga'][i]/Sbase
            barra['P'] = temp3[:]
            barra['Q'] = temp4[:]
        del(temp3)
        del(temp4)
        temp5 = []
        for i in range(nb):
            DadosBarra['tensao'][i] = DadosBarra['tensao'][i]/1000
            DadosBarra['angulo'][i] = DadosBarra['angulo'][i]*(np.pi/180)   # Converte de radianos para graus para radianos?
            DadosBarra['Qshunt'][i] = DadosBarra['Qshunt'][i]/Sbase
            temp5.append(DadosBarra['Qshunt'][i])
            barra['bsh']= temp5[:]
        del(temp5)
        # print(barra['Q'])

        # Convertendo dados da linha

        for i in range(nl):
            if DadosLinhas['tap'][i] == 0:
                DadosLinhas['tap'][i] = 1

        ntn = 0
        ntd = 0
        temp3 = []
        for i in range(nl):
            DadosLinhas['R'][i] = DadosLinhas['R'][i]/100               # Converte de real para pu
            DadosLinhas['Xl'][i] = DadosLinhas['Xl'][i]/100             # Converte de real para pu
            temp3.append(DadosLinhas['Xc'][i]/200)                      # Calculo do Bsh
            DadosLinhas['fi'][i] = DadosLinhas['fi'][i]*(np.pi/180)     # Tranforma o angulo do trafo defasador em rad
            if DadosLinhas['tap'][i] != 1.0:
                ntn += 1                                                # Conta o numero de trafos normais
            elif DadosLinhas['fi'][i] != 0.0:
                ntd += 1                                                # Conta o numero de trafos defasadores

        # --------------------------------------------------------------------------------------------------------------------------------------------------------
        #       Montagem da matriz Y, G e B PARA O FLUXO DE CARGA

        Y = np.zeros((nb,nb),dtype=complex)
        ramo = {'g':None,'b':None,'y':None,'bsh':None}
        temp = []
        temp1 = []
        temp2 = []

        for J in range(nl): 
            temp.append(DadosLinhas['R'][J]/(np.power(DadosLinhas['R'][J], 2)+(np.power(DadosLinhas['Xl'][J], 2))))
            temp1.append(-DadosLinhas['Xl'][J]/(np.power(DadosLinhas['R'][J], 2)+(np.power(DadosLinhas['Xl'][J], 2))))
            temp2.append(temp[J]+temp1[J]*1j)
            ramo['g'] = temp[:]
            ramo['b'] = temp1[:]
            ramo['y'] = temp2[:]
            ramo['bsh'] = temp3[:]
        del(temp)
        del(temp1)
        del(temp2)
        del(temp3)
        #print(ramo['bsh'])
        #print(len(ramo['bsh']))

        for J in range(nl):
            k = int(DadosLinhas['De'][J])
            m = int(DadosLinhas['Para'][J])
            Y[k][m] = -(1/DadosLinhas['tap'][J]) * (np.exp(-DadosLinhas['fi'][J]*1j)) * ramo['y'][J]
            Y[m][k] = -(1/DadosLinhas['tap'][J]) * (np.exp(DadosLinhas['fi'][J]*1j)) * ramo['y'][J]
            Y[k][k] = Y[k][k] + ((1/DadosLinhas['tap'][J])**2)*ramo['y'][J] + 1j*ramo['bsh'][J]
            Y[m][m] = Y[m][m] +                                ramo['y'][J] + 1j*ramo['bsh'][J]

        for i in range(nb):
            Y[i][i] = Y[i][i] + barra['bsh'][i]*1j        #Me atentar se preciso ou nao tirar 1 desse laco

        G = np.real(Y)
        B = np.imag(Y)

        #print(Y)
        #print('--'*30)
        #print(B)

        return[DadosBarra,DadosLinhas,barra,ramo,Sbase,nb,nl,Y,B,G,ntn,ntd]