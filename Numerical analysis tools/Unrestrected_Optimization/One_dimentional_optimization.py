import numpy as np

def fixed_Step_Size(f,x1, p ):

    x2= x1+ p
    f2=f(x2)
    f1=f(x1)
   
    if (f2<f1) :
        while (f2<f1) :
             x1 = x2
             x2 = x2+p
             f1 = f2
             f2 =f(x2)
             
    elif (f2>f1) :
        while (f2<f1) :
            x2 = x1
            x1 =x1-p
            f2 =f1
            f1 = f(x1)
           
    return (x1+x2)/2

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

def exhaustive_method(f ,xs,xf,n):
     
     l=np.linspace(xs,xf,n)
     fonction =[]
     for i in range(0,n):
         fonction.append(f( l[i] ))
     minimum = min(fonction)   
     indice_min = fonction.index(minimum) 
     return l[indice_min]

def dichotomous_method(f,xs,xf,n):
    
    delta=0.001
    l=xf-xs
    x1 =xs+ (l/2) -(delta/2)
    x2 =xs +(l/2) + (delta/2)
    f1=f(x1)
    f2=f(x2)
    for i in range(1,n+1) :
        if f1<f2 :
            xf =x2

        elif f1>f2:
            xs=x1    
        else :
            xs=x1
            xf=x2
        l=xf-xs
        x1 = xs +(l/2) -(delta/2)
        x2 = xs +(l/2) +(delta/2)
        f1=f(x1)
        f2=f(x2)
    return (x1+x2)/2

def intervalhalving_method(f,xs,xf,n):
   import numpy as np

   l=xf-xs
   for i in range(n+1):
       divide =np.linspace(xs,xf,5)
       x1=divide[1]
       x0=divide[2]
       x2=divide[3]
       f1,f0,f2 =f(x1),f(x0),f(x2)
       if f0>f1 and f0<f2 :
           xf = x0
       elif f0<f1 and f0>f2  :
           xs=x0 
       elif f0<f1 and f0<f2 :
           xs=x1
           xf =x2
       l=xf-xs
   return (xs+xf)/2

def fibo(n):
    if n==0 or n== 1 :
        return 1
    else :
        un_2=1
        un_1=1
        for i in range (n,1,(-1)) :
            un=un_1+un_2
            un_2,un_1=un_1,un
    return un  

def fibonacci_method(f,xs,xf,n) :
   
    l0=xf-xs
    l= ( fibo(n-2)/fibo(n) ) * l0
    x1=xs+l
    x2=xf-l 
    f1 ,f2 =f(x1) , f(x2)

    while( x1!=x2 ) :
        if f1<f2 :
            xf =x2
            x2 =x1
            n =n-1
            l0=xf-xs
            l= ( fibo(n-2)/fibo(n) ) * l0
            x1=xs+l
            f1 ,f2 =f(x1) , f1

        elif f1>f2 :
            xs=x1
            x1=x2
            n=n-1
            l0=xf-xs
            l= ( fibo(n-2)/fibo(n) ) * l0
            x2=xf-l
            f1,f2 =f2,f(x2)
    return x1

def goldensection_method(f,xs,xf):
   
    l0 =xf-xs
    
    phi = 0.618
    l= ( 1-phi ) * l0
  
    x1=xs+l
    x2=xf-l 
    f1 ,f2 =f(x1) , f(x2)
   

    while( abs(x2-x1)>0.01  ) :
        if f1<=f2 :
            xf =x2
            x2 =x1
            l0=xf-xs
            l= ( 1-phi ) * l0
            x1=xs+l
            f1 ,f2 =f(x1) , f1
           

        elif f1>=f2 :
            xs=x1
            x1=x2
            l0=xf-xs
            l= ( 1-phi ) * l0
            x2=xf-l
            f1,f2 =f2,f(x2)
           
    return x2

