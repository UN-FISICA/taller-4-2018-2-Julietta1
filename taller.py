# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 19:15:25 2018

@author: salvi
"""
import numpy as np
from math import pi
from scipy import optimize
class Derivada:
    def __init__(self,f,metodo="adelante",dx=0.001):
        self.f=f
        self.metodo=metodo
        self.dx=dx
    def calc(self,x):
        if self.metodo=="adelante":
            der=(self.f(x+self.dx)-self.f(x))/self.dx
            return(der)
        if self.metodo=="central":
            der=((self.f(x+(self.dx)/2)-self.f(x-(self.dx)/2))/self.dx)
            return (der)
        if self.metodo=="extrapolada":
            derf1=(self.f(x+(self.dx)/2)-self.f(x-(self.dx)/2))/self.dx
            derf2=(self.f(x+(self.dx)/4)-self.f(x-(self.dx)/4))/(self.dx/2)
            return (4*derf2-derf1)/3
        if self.metodo=="segunda":
            der2=((self.f(x+self.dx)+self.f(x-self.dx)-2*self.f(x))/(self.dx)**2)
            return(der2)
        
            
            
class Zeros:
    def __init__(self,f,metodo,error=1e-4,max_iter=100):
        self.f=f
        self.metodo=metodo
        self.error=error
        self.max_iter=max_iter
    def zero(self,vi):
        if self.metodo=="newton":
            derf=Derivada(self.f)
            x=float(vi)
            for i in range(self.max_iter):
                xi=x-(self.f(x)/derf.calc(x))
                x=xi
                if self.f(x)==0:
                    return(xi)
            if self.f(xi)<=self.error:

                return(xi)
            else:
                return("No se encontró una raiz, pruebe con un número mayor de iteraciones o cambiando el valor inicial")
        if self.metodo=="bisectriz":
            a=vi[0]
            b=vi[1]
            for i in range(self.max_iter):
                xi=(a+b)/2
                if abs(self.f(b))<abs(self.f(a)):
                    a=xi
                else:
                    b=xi
                if self.f(xi)==0:
                    return(xi)
            if abs(self.f(xi))<=self.error:
                return(xi)
            else:
                return("No se encontró una raíz, pruebe con un número mayor de iteraciones o cambiando el valor inicial")
                
        if self.metodo=="interpolacion":
            a=vi[0]
            b=vi[1]
            y0=self.f(a)
            y1=self.f(b)
            for i in range(self.max_iter):
                c=a-((y0*(b-a))/(y1-y0))
                if self.f(c)==0:
                    return(c)
                if abs(self.f(a))<abs(self.f(b)):
                    if (self.f(c)<0 and self.f(b)<0) or (self.f(c)>0 and self.f(b)>0):
                        b=c
                    else:
                        a=c
                    
                else:
                    if (self.f(c)<0 and self.f(a)<0) or (self.f(c)>0 and self.f(a)>0):
                        a=c
                    else:
                        b=c
            if abs(self.f(c))<=self.error:
                return(c)
            else:
                 return("No se encontró una raíz, pruebe con un número mayor de iteraciones o cambiando el valor inicial")    
        if self.metodo=="newton-sp":
            raiz=optimize.newton(self.f,vi,maxiter=self.max_iter,tol=self.error)
            return(raiz)
        if self.metodo=="fsolve-sp":
            raiz=optimize.fsolve(self.f,vi,xtol=self.error)
            return(raiz)
        if self.metodo=="brentq-sp":
            raiz=optimize.brentq(self.f,vi[0],vi[1],maxiter=self.max_iter,xtol=self.error)
            return(raiz)
            
if __name__=="__main__":
    a=Derivada(np.sin).calc(pi)
    b=Derivada(np.sin,"central").calc(pi)
    c=Derivada(np.sin,"extrapolada",0.00001).calc(pi)
    d=Derivada(np.cos,"segunda").calc(pi)
    
    e=Zeros(np.sin,"newton").zero(3)
    f=Zeros(np.cos,"bisectriz",1e-6).zero((1.0,2.0))
    g=Zeros(np.log,"interpolacion",1e-6,90).zero((0.5,1.5))
    h=Zeros(np.sin,"newton-sp").zero(3.3)
    i=Zeros(np.cos,"fsolve-sp").zero(1)
    j=Zeros(np.log,"brentq-sp").zero((0.5,1.5))
    
    print("f'(pi)=",a)#imprime la derivada de seno evaluada en pi con el método adelante por defecto y el dx=0.001 por defecto
    print("f'(pi)=",b)#imprime la derivada de seno evaluada en pi con el método central y el dx=0.001 por defecto
    print("f'(pi)=",c)#imprime la derivada de seno evaluada en pi con el metodo de la derivada extrapolada y dx=0.00001
    print("f''(pi)=",d)#imprime la segunda derivada de coseno evaluada en pi
    print("Raíz numérica encontrada en x=",e)#imprime la raiz de seno con valor inicial 3 por el método de newton, un error máximo de 1e-4 y 100 iteraciones máximas por defecto
    print("Raíz numérica encontrada en x=",f)#imprime la raiz de coseno en el intervalo (1,2) por el metodo de bisectriz con un error máximo de 1e-6 
    print("Raíz numérica encontrada en x=",g)#imprime la raiz de ln en el intervalo (0.5,1.5) por el método de interpolacion con un error máximo de 1e-6 y 90 iteraciones máximas
    print("Raíz numérica encontrada en x=",h)#imprime la raiz de seno usando la función de newton definida en scipy
    print("Raíz numérica encontrada en x=",i)#imprime la raiz de coseno usando la función fsolve definida en scipy
    print("Raíz numérica encontrada en x=",j)#imprime la raiz de ln usando la función brentq definida en scipy