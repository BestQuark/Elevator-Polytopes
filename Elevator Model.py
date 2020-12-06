#!/usr/bin/env python
# coding: utf-8

# # Optimización lineal sobre el politopo de Birkhoff aplicado a modelo sencillo de elevadores

# ### En la siguiente celda se define el modelo de elevadores planteado en [1] y se definen un par de funciones para correr el programa

# In[171]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import linprog
import warnings

    
# Esta función dibuja el elevador

def printElevador(numero_Personas, numero_Elevadores, numero_Pisos, posicion_Personas, posicion_Elevadores):
    
    print("\n \n \n Modelo hiperealista del Burj Khalifa \n")
    
    
    n_personas_repetidas = np.amax(np.array([np.count_nonzero(posicion_Personas == i) for i in posicion_Personas]))   
    n_elevadores_repetidos = np.amax(np.array([np.count_nonzero(posicion_Elevadores == i) for i in posicion_Elevadores]))       
    for i in reversed(range(1,numero_Pisos+1)):
        
        k = len(str(i))
        n_perPiso = np.count_nonzero(posicion_Personas == i)
        elevs = ''
        for j in posicion_Elevadores:
            if j==i:
                elevs = elevs+'|████|'
            else:
                elevs = elevs+'|    |'
        print( "|"+'_'*5 +'_'*4*n_personas_repetidas + '_'*numero_Elevadores*6 + "| \n")
        print( "|"+' '*(4*(n_personas_repetidas-1))+"Piso "+str(i)+' '*(3-k) + "|"  + elevs +"| \n")
        print( "|"+' '*(4+4*(n_personas_repetidas-n_perPiso))+' ●  '*n_perPiso + "|" +elevs+"| \n")
        print( "|"+' '*(4+4*(n_personas_repetidas-n_perPiso))+'¬█¬ '*n_perPiso + "|" +elevs+"| \n")
        print( "|"+' '*(4+4*(n_personas_repetidas-n_perPiso))+' n  '*n_perPiso+ "|" +elevs+"| \n")
    print("|"+'_'*(5+4*n_personas_repetidas+6*numero_Elevadores)+"|")
    
# Esta función resuelve el problema LP

def solElevador(numero_Personas, numero_Elevadores, numero_Pisos, posicion_Personas, posicion_Elevadores, metodo):
    warnings.filterwarnings('ignore')
    if numero_Personas>numero_Elevadores:
        numero_Esperando = numero_Personas-numero_Elevadores
        posicion_Esperando = posicion_Personas[numero_Esperando:]
        print("Hay ", numero_Esperando , " persona(s) esperando a que se desocupen los ascensores.")
        numero_Personas = numero_Elevadores
        posicion_Personas = posicion_Personas[:numero_Elevadores]
       
    l = numero_Elevadores**2
    c = np.zeros(l)
    for i,posP in enumerate(posicion_Personas):
        c[numero_Elevadores*i:numero_Elevadores*(i+1)] = np.abs(posicion_Elevadores - posP)   
    
    A = np.zeros([2*numero_Elevadores,l])
    for i in range(numero_Elevadores):
        A[i, numero_Elevadores*i:numero_Elevadores*(i+1)] = np.ones(numero_Elevadores)
        A[i+numero_Elevadores, i::numero_Elevadores] = np.ones(numero_Elevadores)
    
    b = np.ones(2*numero_Elevadores)
    return linprog(c, A_eq=A, b_eq=b, bounds=(0,1), method=metodo)

def simuladorElevador():
    numero_Elevadores = int(input("Ingrese el número de elevadores: "))
    numero_Personas = int(input("Ingrese el número de personas: "))
    numero_Pisos = int(input("Ingrese el número de pisos del edificio: "))
    print("\n"+"Métodos: \n")
    print("1: Simplex Method")
    print("2: Interior Point Method")
    print("3: Revised Simplex Method")
    n_metodo = int(input("Ingrese el método que desee usar:"))
    metodo = 0
    if n_metodo == 1:
        metodo = 'simplex'
    elif n_metodo == 2:
        metodo = 'interior-point'
    elif n_metodo == 3:
        metodo = 'revised simplex'
    if (metodo==0) or (numero_Personas<=0) or (numero_Pisos <= 0) or (numero_Elevadores <= 0):
        print("\n"+"Los valores ingresados son inválidos")
    else:
        posicion_Personas = np.array([np.random.randint(1,numero_Pisos+1) for i in range(numero_Personas)])
        posicion_Elevadores = np.array([np.random.randint(1,numero_Pisos+1) for i in range(numero_Elevadores)])        
        printElevador(numero_Personas, numero_Elevadores, numero_Pisos, posicion_Personas, posicion_Elevadores)           
        linPsolution = solElevador(numero_Personas, numero_Elevadores, numero_Pisos, posicion_Personas, posicion_Elevadores,metodo)
        M = np.reshape(linPsolution.x, [numero_Elevadores,numero_Elevadores])
        personas_orden = posicion_Personas[:numero_Elevadores]
        elevadores_orden = np.matmul(M,np.array(range(1,numero_Elevadores+1)).transpose())
        for i in range(min(numero_Elevadores,numero_Personas)):
            print("El elevador "+ str(int(elevadores_orden[i]))+ " recoge a la persona del piso "+ str(personas_orden[i]))

        print("\n"+"Datos usados: \n")
        print("Posición de las personas recogidas: \n",personas_orden, "\n")
        print("Posición de todas las personas: \n",posicion_Personas, "\n")
        print("Posición de los elevadores: \n ", posicion_Elevadores, "\n")
        print("Matriz de permutación de los elevadores: \n", M.astype(int))


# In[172]:


simuladorElevador()


# ## Referencias
# 
# [1] Mehlum, M. (2012). Doubly stochastic matrices and the assignment problem (Master's thesis).
