from math import* 
import numdifftools as nd 
import numpy as np

#Gradient Descent Method 
def norme (x):
  return sqrt(np.sum(x**2))

def accelerated_sep_size(f,x1 ,epsilon):
    
    h= epsilon/2
    x2 = x1+h
    f2=f(x2)
    f1=f(x1)
    if f2>f1 :
        pas=h
        while f2>f1 :
             pas = pas*2
             x2 =x1 
             x1 =x1-pas
             f2 = f1
             f1 = f(x1)
        while abs(x2-x1) >= epsilon :
            pas = h 
            while (f1>f2) and abs(x2-x1)>=epsilon and x1<x2 :
                pas =pas*2
                x1=x1+pas
                f1=f(x1)
            if abs(x2-x1)<=epsilon :
                return (x1+x2)/2
            elif x1>x2 :
                f1,f2 =f2,f1
                x1,x2 =x2,x1

            pas= h

            while (f1<f2) and abs(x2-x1)>=epsilon and x1<x2 :
                pas=pas*2
                x2=x2-pas
                f2=f(x2)
            if abs(x2-x1)<=epsilon :
                return(x1+x2)/2 
            elif x1>x2 :
                f1,f2=f2,f1 
                x1,x2=x2,x1     
    elif f2 <f1 :
        pas=h
        while f2<f1 :
            pas = pas*2
            x1 =x2
            x2=x2+pas  
            f1 =f2
            f2 =f(x2)
        
        while abs(x2-x1)>=epsilon :
         pas=h
         while (f2>f1) and (abs(x2-x1)>= epsilon) and x2> x1 :
                pas*=2
                x2=x2-pas
         if (abs(x2-x1)<= epsilon):
            return (x1+x2)/2
         elif x2<x1 :
            f1,f2=f2,f1
            x2,x1=x1,x2
            
         while (f2<f1) and abs(x2-x1)>=epsilon and x1<x2 :
                pas=pas*2
                x1=x1+pas
                f1=f(x1)
         if abs(x2-x1)<=epsilon :
            return (x1+x2)/2
         elif x1>x2 :
            f2,f1 =f1 ,f2 
            x1 ,x2=x2,x1 

    elif f2 ==f1:
        return (x1+x2)/2

def Gradient_descent (f,x,epsilon):
  path =[x]
  df=nd.Gradient(f)
  while (norme(df(x)) > epsilon) :
    g = lambda alpha : f(x- alpha *df(x)) 
    alpha=accelerated_sep_size(g, 1 ,epsilon)
    x = x-alpha* df(x)
    path+=[x]
    
  return x , path

#Conjugate Gradient Method 
def Conjugate_Gradient(f,x ):
    
    path=[x]
    n=len(x)

    df=nd.Gradient(f)
    grd=df(x)
    d=-grd

    Hessfunc=nd.Hessian(f)
    Q=Hessfunc(np.zeros(n))
    
    for k in range(0,n):
        
        alpha = ( (d.T) @ d) / ( (d.T)@ Q @ d )
        
        x = x + alpha*d
        path+=[x]
        
        grd=df(x)
        beta= (grd.T)@ Q@d / (d.T)@ Q @ d
        
        d = -grd +beta * d
        
    return x ,path

#Newton METHOD    
def definie_positive(A) :
 Spectre=np.linalg.eigvals(A)
 for e in (Spectre):
     if e<=0 :
         return False
 return True

def symetric(A):
    n=len(A)
    T_A=np.transpose(A)
    
    for i in range(n):
        for j in range(n):
            if A[i][j]!=T_A[i][j]:
               return False
    
    return True

def Newton (f,x,epsilon ):
   
   path =[x]
   n=len(x)
 
   df=nd.Gradient(f)
   grd=df(x)
  
   Hessfunc=nd.Hessian(f) 
   H=Hessfunc(x)
   
   d=- (np.linalg.inv(H))@grd
   
   k=0
   while (norme(d) > epsilon) :
     
     g = lambda alpha : f(x + alpha *d ) 
     alpha=accelerated_sep_size(g, 1 ,epsilon)  
     x = x+ alpha*d
     path+=[x]
     
     H=Hessfunc(x)
     grd= df(x)

     if (symetric(H) and  definie_positive(H) ):
         d=-(np.linalg.inv(H))@grd
     else :
         eps=min(Spectre=np.linalg.eigvals(H))
         I=np.identity(n)

         d= - np.linalg.inv(eps*I+H) @grd
     
     
     k+=1

   return x,path

#Quasi Newton METHOD
def armijo(phi,alpha0,nu,epsilon,d,grd):
    
    phi_approx = lambda alpha : phi(0)+ epsilon* (d.T @ grd) * alpha

    alpha=alpha0 
    
    if  phi( alpha )<phi_approx(alpha) :
        #forward research = back tracking
      while phi( alpha )<phi_approx(alpha) :
            ancienne_valeur =alpha
            alpha=nu*alpha
      return ancienne_valeur 
    
    
    else :
        #backward research = back tracking
       while phi( alpha ) >= phi_approx(alpha) :
            ancienne_valeur =alpha
            alpha= alpha/nu            
       return ancienne_valeur

def Quasi_Newton(f,x,tolerance):
 path =[x]
 #les parametres d'armijo 
 epsilon=tolerance
 nu=2
 alpha = epsilon 
 
 n=len(x)
 df =nd.Gradient(f)
 grd=df(x)
 H=np.eye(n,n)

 while ( norme(grd) > tolerance ):
     d = - H@grd
     phi= lambda alpha : f(x+alpha*d)  

     
     alpha_min  =armijo(phi,alpha,nu,epsilon,d,grd)
     
     x=x+alpha_min*d
     path+=[x]

     y=-grd
     grd=df(x)  
     y+=grd
     A= (alpha_min* (d @ d.T) )/ (d.T @ y )
     B= - ( H @ y ) @ (( H@y ).T) / ( y.T @ H @ y )
     H= H +A + B
     
 return x,path
