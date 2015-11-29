import matplotlib.pyplot as plt
import numpy as np

def f_oben(x):
    return 1.5*(2*x-1)**2+0.25

def f_unten(x):
    return -1.5*(2*x-1)**2-0.25

thirdpoint1 = [0.85, 0.5]
thirdpoint2 = [0.15, 0.5]

A = np.array([[1,0, 0],[1, 1, 1], [1, thirdpoint1[1], thirdpoint1[1]**2]])
y = np.array([1,1,thirdpoint1[0]])
c1 = np.linalg.solve(A,y)

A = np.array([[1,0, 0],[1, 1, 1], [1, thirdpoint2[1], thirdpoint2[1]**2]])
y = np.array([0,0,thirdpoint2[0]])
c2 = np.linalg.solve(A,y)

x = np.arange(0, 1.05, 0.05)



def f_rechts(y):
    return c1[0]+c1[1]*y+c1[2]*y**2

def f_links(y):
    return c2[0]+c2[1]*y+c2[2]*y**2

def boundary(chi, nu):
    if type(nu) is not np.ndarray and type(chi) is not np.ndarray and nu == 0 and chi == 0:
        return 0, f_unten(0)
    if type(nu) is not np.ndarray and nu == 0:
        return chi, f_unten(chi)
    if type(nu) is not np.ndarray and nu == 1:
        return chi, f_oben(chi)
    if chi == 0:
        return f_links(nu), nu*f_oben(0)+(1-nu)*f_unten(0)
    if chi == 1:
        return f_rechts(nu), nu*f_oben(1)+(1-nu)*f_unten(1)

def grid(chi, nu):
    chiarr = np.array([1-chi, chi])
    nuarr = np.array([1-nu, nu])
    rnuarrx = np.array([boundary(0,nu)[0], boundary(1, nu)[0]])
    rchiarrx = np.array([boundary(chi,0)[0], boundary(chi,1)[0]])
    rmatrx = np.array([[boundary(0,0)[0], boundary(0,1)[0]],[boundary(1,0)[0],boundary(1,1)[0]]])
    
    rnuarry = np.array([boundary(0,nu)[1], boundary(1, nu)[1]])
    rchiarry = np.array([boundary(chi,0)[1], boundary(chi,1)[1]])
    rmatry = np.array([[boundary(0,0)[1], boundary(0,1)[1]],[boundary(1,0)[1],boundary(1,1)[1]]])

    x = chiarr.dot(rnuarrx)+rchiarrx.dot(nuarr)-chiarr.dot((rmatrx.dot(nuarr)))

    y= chiarr.dot(rnuarry)+rchiarry.dot(nuarr)-chiarr.dot((rmatry.dot(nuarr)))

    return x, y


def first2():
    plt.plot(boundary(x, 0)[0], boundary(x,0)[1], color="black")
    plt.plot(boundary(x,1)[0], boundary(x,1)[1], color="black")

    plt.plot(boundary(0, x)[0], boundary(0,x)[1], color="black")
    plt.plot(boundary(1,x)[0], boundary(1,x)[1], color="black")
    plt.show()

def x_interp2():
    plt.plot(boundary(x, 0)[0], boundary(x,0)[1], color="black")
    plt.plot(boundary(x,1)[0], boundary(x,1)[1], color="black")


    xbegin, ybegin = np.array(boundary(x,0)[0]),np.array(boundary(x,0)[1])
    xend, yend = np.array(boundary(x,1)[0]),np.array(boundary(x,1)[1])

    for r in x:
        xnow = r*xbegin+(1-r)*xend
        ynow = r*ybegin+(1-r)*yend
        plt.plot(xnow, ynow, color="red")

    for r in x:
        plt.plot([r, r], [f_unten(r), f_oben(r)], color="black")

    plt.show()

def y_interp2():

    plt.plot(boundary(0, x)[0], boundary(0,x)[1], color="black")
    plt.plot(boundary(1,x)[0], boundary(1,x)[1], color="black")

    for i in range(0,len(x)):
        a, b = boundary(0,x)[0][i],boundary(0,x)[1][i]
        c, d = boundary(1,x)[0][i],boundary(1,x)[1][i]
        plt.plot([a,c], [b,d], color="black")

    xbegin, ybegin = np.array(boundary(0,x)[0]),np.array(boundary(0,x)[1])
    xend, yend = np.array(boundary(1,x)[0]),np.array(boundary(1,x)[1])

    for r in x:
        xnow = r*xbegin+(1-r)*xend
        ynow = r*ybegin+(1-r)*yend
        plt.plot(xnow, ynow, color="red")

    plt.show()

def whole2():
    plt.plot(boundary(x, 0)[0], boundary(x,0)[1], color="black")
    plt.plot(boundary(x,1)[0], boundary(x,1)[1], color="black")
    plt.plot(boundary(0, x)[0], boundary(0,x)[1], color="black")
    plt.plot(boundary(1,x)[0], boundary(1,x)[1], color="black")

    for r in x:
        xnow = [grid(r, t)[0] for t in x]
        ynow = [grid(r, t)[1] for t in x]
        plt.plot(xnow, ynow, color="black")
        xnow = [grid(t, r)[0] for t in x]
        ynow = [grid(t, r)[1] for t in x]
        plt.plot(xnow, ynow, color="black")

    plt.show()



first2()
x_interp2()
y_interp2()
whole2()