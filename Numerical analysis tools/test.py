import numpy as np 
import matplotlib.pyplot as plt


#one dimensional Optimization
import Unrestrected_Optimization.One_dimentional_optimization
def g(x):
    return x*(x-1.5)
xs=0.0
xf=2.0
n =6
min=Unrestrected_Optimization.One_dimentional_optimization.fibonacci_method(g,xs,xf,n)

x=np.linspace( 0 , 2, 30 )
fig, ax = plt.subplots()
ax.plot(x, g(x), label='g(x)')
ax.grid(True)
ax.set_xlabel('axe des x')
ax.set_ylabel('axe des y')
ax.set_title("One dimensional Optimazation via unrestrected methods")
ax.plot(min,g(min),'r*', label="min")
ax.legend()
ax.annotate ( f" minimum local = ( {min} , {g(min)} ) ", xy=(min, g(min)), xytext=(0.5, 0.01), arrowprops=dict(facecolor='black', shrink=0.05) )
plt.show()
plt.close()

#Multidimensional Optimization
import Unrestrected_Optimization.Multivariable_optimization
f= lambda x:   0.5* ( x[0]**2+x[1]**2 ) 
x=np.array([2,1])
epsilon =1e-15
min,path=Unrestrected_Optimization.Multivariable_optimization.Gradient_descent (f,x,epsilon)


fig,ax =plt.subplots(subplot_kw={"projection": "3d"} )
x_axe=np.linspace(-3,3,100)
y_axe=np.linspace(-3,3,100)
X,Y=np.meshgrid(x_axe,y_axe)
Z=np.array([X,Y])
F=f(Z)

ax.plot_surface(X, Y, F ,antialiased=True,cmap= 'inferno')
ax.scatter( min , min,  f( np.array( [min,min] ) ) ,color="red",s=200)
ax.view_init(azim=30, elev= 15)
ax.set_title("Multidimensional Optimazation via unrestrected methods")
plt.show()
plt.close()

plt.contour(X,Y,F,10)
plt.title("Contour Plot")
plt.show()
plt.close()



#Comparison between methods 
  
  #Create 4 subplots 
fig,axs=plt.subplots(2,2,figsize=(10,10))
fig.suptitle("Comparison between Multidimensional methods" )

  #For each subplot we create title and contour plot
Methods=["Gradient Descent","Conjugate Gradient ","Newton","Quasi Newton "]  
i=0
for row in range(2):
    for col in range(2):
        #axs[row,col].annotate(f'{Methods[i]} Method',(0.5, 0.5), transform=axs[row, col].transAxes)
        axs[row,col].set_title(f' {Methods[i]} Method ')
        axs[row,col].contour(X,Y,F,10)
        i+=1



  #Gradient Descent
min0,path0=Unrestrected_Optimization.Multivariable_optimization.Gradient_descent (f,x,epsilon)
print(path0)
         #make "patth0" an np.array type
path0=np.array(path0)
        #path0[:,0] _example :
""" L=np.array([[1,2],[7,8]])             L=[[1 2]
                                           [7 8]]
    print(L)                      >>>>>  L[:,0] =[1 7] 
    print(L[:,0])
"""
   
F0=[]
for e in path0 :
    F0+=[f(e)]
    
axs[0,0].scatter(path0[:,0], F0 ,c='black' ,marker='o')
axs[0,0].plot( path0[:,0] , F0 , color = 'red', linestyle = 'solid') 
  
  #Conjugate Gradient
min1,path1=Unrestrected_Optimization.Multivariable_optimization.Conjugate_Gradient(f,x )
path1=np.array(path1)
F1=[]
for e in path1:
    F1+=[f(e)]
axs[0,1].scatter( path1[:,0] , F1, c='black' ,marker ='o' )
axs[0,1].plot( path1[:,0] ,F1 , color='red' ,linestyle ='solid' )
  
   #Newton
min2,path2=Unrestrected_Optimization.Multivariable_optimization.Newton (f,x,epsilon )
path2=np.array(path2 )
F2=[]
for e in path2 :
    F2+=[f(e)]
axs[1,0].scatter(path2[:,0] , F2 ,c='black',marker='o' )
axs[1,0].plot(path2[:,0], F2, color='red',linestyle='solid')
   #Quasi Newton
min3,path3=Unrestrected_Optimization.Multivariable_optimization.Quasi_Newton(f,x,epsilon)
path3=np.array(path3)
F3=[]
for e in path3 :
    F3+=[f(e)]
axs[1,1].scatter(path3[:,0], F3 , c='black',marker='o' )
axs[1,1].plot(path3[:,0] , F3 ,color='red' , linestyle='solid')

plt.show()
plt.close()


