import numpy as np


def testPrecision():
    print(help(np.finfo))
    
    for dtype in [float, np.double, np.single, np.half]:
        print(dtype.__name__, np.finfo(dtype))


def equalityTestError():
    a = np.sqrt(2.0)
    print(f"a={a}, type(a)={type(a)}")
    
    b = a*a
    print(f"b={b}")
    
    print(b == 2.0)
    print(b-2.0)
   
def equalityTestErrorFixed():
    a = np.sqrt(2.0)
    print(f"a={a}, type(a)={type(a)}")
    b = a*a
    print(f"b={b}")
    print(my_isclose(a*a, b))
    print(np.isclose(a*a, b))

def my_isclose(x, y, tol=1.0e-9):
    return abs(x - y) < tol


def f(x):
    return x*(np.sqrt(x+1)-np.sqrt(x))

def g(x):
    return x/(np.sqrt(x+1)+np.sqrt(x))



def mainTestforfandgfunctions():
    x = np.float64(500.0)
    y = np.float32(500.0)
    z = np.float16(500.0)
    print("float 64 prevision f(x):", f(x))
    print("float 64 prevision g(x):", g(x))
    print("difference: ", abs(f(x)-g(x)))
    print("\t---")
    print("float 32 prevision f(x):", f(y))
    print("float 32 prevision g(x):", g(y))
    print("difference: ", abs(f(y)-g(y)))
    print("\t---")
    print("float 16 prevision f(x):", f(z))
    print("float 16 prevision g(x):", g(z))
    print("difference: ", abs(f(z)-g(z)))
    

def piSquared(n):
    piS = 0.0
    k = 1
    while k < (0.5+n):
        piS += (k)**(-2)
        k+=1
    return 6.0*piS

def piSquared_v2(n): #this doesn't work for some reason.
    output=np.float64(0)
    for i in range(n,1, -1):
        output += i**(-2)
    return 6.0*output

def alternative(n): #better approach
    piS = 0.0
    k = n
    while k > (0.5):
        piS += (k)**(-2)
        k-=1
    return 6.0*piS

def piSquaredTest():
    piSquare = np.float64(np.pi**2)
    n = [1.e6, 1.e7, 1.e8, 1.e9]
    for i in n:
        x = np.float64(piSquared(i))
        print("For increment function: ",i, "\tResult: ",x, "\tDifference: ", abs(piSquare-x))
        
        #y = np.float64(piSquared_v2(i))
        y = np.float64(alternative(i))
        print("For decrement function: ",i, "\tResult: ",y, "\tDifference: ", abs(piSquare-y))
        
   

#piSquaredTest()

print(np.finfo(np.float32))

def largest(b,t,l,u): #write code to compute this
    y = np.float64(0) #output
    b = np.float64(b)
    t = int(t)
    l = np.float64(l) 
    u = np.float64(u)     
    
    #TODO: compute y=0.11..1 x 2^u
    
    for i in range(t):
        y += b**(-i-1) # calculates the value for each binary -> calculates the fraction part.
    return b**u*y #2^126 * digits
    
    

def mainLargest():
    print(f"The largest number is {largest(2,24,-126,128)}.")

mainLargest()

