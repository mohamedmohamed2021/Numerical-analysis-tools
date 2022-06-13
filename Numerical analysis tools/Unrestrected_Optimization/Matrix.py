from fractions import Fraction
import numpy as np

# GAUSSIAN RESOLUTION
def eliminer(j,A,pivot):
    for l in range(j+1,len(A)) :
        p=A[l][j]
        for c in range(j,len(A[0])):
            A[l][c] = A[l][c] - Fraction( (p*A[j][c]) / pivot )
            
def triangulariser (A):
    for j in range (len(A)-1) : 
       eliminer(j,A,A[j][j])

def backward(A) :
    n=len(A)
    X=np.zeros(n)
    X[(n-1)] = A[(n-1)][n] / A[(n-1)][(n-1)]
    for i in range(n-2,-1,-1):
       p=A[i][n]
       S= 0
       
       for k in range (n-1,i,-1) :
           S+=A[i][k]*X[k]
      
       X[i] =  (p-S)/A[i][i]

    return X   

def Gauss_system_resolution(A,B):
    n=len(A)
    """we will build the augmented matrix"""
    B=np.reshape(B,(n,1))
    A=np.append(A, B,axis= 1) 
    
    
    
    triangulariser(A)
   
    X=backward(A)
    return X

#GAUSSIAN METHOD TO INVERSE MATRIX
def elimination_inverse(j,A,pivot):
     for l in range(j-1,-1,-1):
       
        p=A[l][j]
       
        for c in range(j,len(A[0])):
            
            A[l][c] = A[l][c] - Fraction( (p*A[j][c]) / pivot )
           
def triangularisation_inverse(A):
    for j in range(len(A)-1,0,-1):
        elimination_inverse(j,A,A[j][j])

def inverse (A) :
  n=len(A)
  I=np.identity(n)
  """we will build an augmented A matrix concatenated with Identity matrix"""
  A=np.concatenate((A,I),axis=1)

 
  triangulariser (A)
 
  triangularisation_inverse(A)

  for i in range(0,n):
      p=A[i][i]
      for j in range(i,len(A[0])):
           A[i][j]= Fraction(A[i][j]/p)
  return A[:,n:2*n]

#LU DECOMPOSITION
def LU_decomposition(A):
    n=len(A)
    L=np.identity(n)
    U =np.zeros((n,n))
    if A[0][0]== 0 :
       PIVOT = pivot (A,0)[0]
       LIGNE_PIVOT = pivot (A,0)[1]
       if PIVOT !=0 :
           echange_ligne(LIGNE_PIVOT,0,A)
    else :
        for j in range(0,n-1):
            for l in range(j+1,n):
                p = A[l][j]
                A[l] =  A[l]- Fraction(p/A[j][j])*A[j]
                L[l][j] =Fraction(p/A[j][j])
                
    U=A
    return L,U

#Choleski's Decomposition 
def Cholesky_decomposition(A):
   n= len(A)
   L=np.zeros((n,n))
   #go all over rows 
   for l in range(n):
       #go all over columns
       for c in range(l,n):
            if (c == l ):
                 
              if (l== 0) :
                     L[0,0]=np.sqrt(A[0,0])
                     
              else:
                 s=0
                 for k in range(0,l):
                   s+= ( L[l,k] ) **2 

                
                 L[l,l]= np.sqrt(A[l,l]-s)
                
            
            else:
              if l==0 :
                    L[c,0] =A[0,c]/L[0,0]
              else:
                s=0
                for k in range(0,l):
                    s+=L[l,k]*L[c,k]
                
                
                
                L[c,l]=( A[l,c] -s) / (L[l,l]) 
               
    
   return L,np.transpose(L)

# LU SYSTEM RESOLUTION

"""verifier si la matrice est définie positive"""
def définie_positive(A) :
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

def forward_substitution(L,B):
    n=len(B)
    Y=np.zeros(n)
    Y[0]=B[0]/L[0,0]
    for k in range(1,n):
        s=0
        for i in range(0,k):
            s+=L[k,i] * Y[i]
        
        Y[k]=(B[k]-s)/L[k,k] 
        
    return Y

def backward_substitution(U,Y):
    n=len(U)
    X=np.zeros(n)
    X[n-1]=Y[n-1]/U[n-1][n-1]
    for k in range(n-2,-1,-1):
        s=0
        for i in range( (k+1) ,n):
            s+= U[k,i] * X[i]
        X[k] = (Y[k]-s) / U[k,k]
    return X        

def LU_system_resolution(A,B):
    if symetric(A) and définie_positive(A)  :
        L,U = Cholesky_decomposition(A)
    else :
        L,U = LU_decomposition(A)    

    Y= forward_substitution(L,B)
    X= backward_substitution(U,Y)
    return X
