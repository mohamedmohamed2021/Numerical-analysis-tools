import math

#POUR NEWTHON on a besoin de l'expression de la dérivée et la dérivée seconde
def Newton_method(f,fp,fpp,xs) :

   xk = xs
   epsilon=1e-5
   while abs(fp(xk))>=epsilon :
     xk = xk -fp(xk)/fpp(xk)
   return xk 

def Quasi_Newton_Method (xs,f):
   h=1e-2
   xk = xs 
   epsilon =1e-3
   fg , fd = f(xk-h) , f(xk+h)
   fc= f(xk)
   
   while abs(( fd - fg ) / (2*h))>= epsilon :
     xk = xk - ( h* ( fd - fg ) ) / (2*( fd + fg - 2* fc ) )
     fg , fd = f(xk-h) , f(xk+h)
     fc= f(xk)
     
   return xk

#POUR Secante_Method on a besoin de l'expression de la dérivée
def Secante_Method (f,fp,xs):
    """ xs: l'abscice le plus proche du minimume depuis la gauche"""
    A = xs
    epsilon=1e-3
    h=1e-2
    while ( fp(h)<0 ):
        A = h
        h*=2
    B=h
    fA , fB=fp(A) , fp(B)
    xk =  A- ( fA * (B-A) ) / (fB-fA)
    fk =fp(xk)
    while abs( fk ) > epsilon :  
        if fk<0 :
            A=xk
            fA=fp(A)
            
        else :
            B=xk
            fB=fp(B)   
        xk= A- ( fA * (B-A) ) / (fB-fA)
        fk =fp(xk) 
    return xk   

