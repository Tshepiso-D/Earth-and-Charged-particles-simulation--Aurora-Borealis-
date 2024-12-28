from vpython import*
## Motion of a charged particle in the (dipole) magnetic field of a circular coil of current
## Based on the original code from Roger Fearick, September 2008
## Updated to work with Python 3 and VPython 7

## Define the circular current-carring loop
def makeloop(rloop,loopradius):
    dtheta=pi/100       # loop is made of segments of this size
    loop = curve(color = vector(1.,.7,.2), radius = 0.01) # visualize loop
    for theta in arange(0.0, 2.0 * pi, dtheta):
        loop.append(rloop + vector(loopradius*sin(theta), 0.0, loopradius*cos(theta)))
    loop.append(loop.point(0)['pos']) # close loop
    return loop

loop = makeloop(vector(0.0,0.0,0.0), 0.3)

## Calculate the magnetic field at r due to the loop
## Note: (mu0/4pi)I arbitrarily set to 10.0
def Field(r, loop):
    B = vector(0,0,0);
    loop_positions = loop.slice(0,-1)
    # print(loop_positions)
    for i in range(len(loop_positions)-1):
        dR = (loop_positions[i]['pos']+loop_positions[i+1]['pos'])/2 - r     # distance r to loop segment
        dI = loop_positions[i]['pos'] - loop_positions[i-1]['pos']           # delta I
        B += cross(dI, dR) / mag(dR)**3
    return 10*B

## Display the magnetic field at a few positions (on the x-y plane)
showfields=1
if showfields:
    #print Field(r0)
    for y in arange(-1.3,1.4,0.2):
        for x in arange(-1.3,1.4,0.2):
            rx=vector(x,y,0)
            B=Field(rx,loop)
            B=B/100.0
            c=color.red
            if mag(B)>0.3:
                B*=0.3/mag(B)
                c=color.magenta
            arrow(pos=rx,axis=B,color=c)

## ============================================

## Model the motion of a charged particle in the magnetic field
            
sphere(radius=0.6,color=color.cyan)
## constants
q = 10.0
m = 1.0
## initial conditions
t = 0.0

r=vector(0.0,4.0,0.0)
v=vector(0.0,-0.5,0.0)
trail=curve(pos=r, radius=0.005, color=color.green)
KE = gcurve()
## integrate
dt = 0.0005
while t < 30:
    # update clock time
  t= t + dt
 # Bfield at r
  B = Field(r,loop)
 # update motion
  magB = mag(B)
  Bhat = B/magB
  f = q*magB/m #effective cyclotron frequency
  theta = f*dt
  sintheta = sin(theta)
  costheta = cos(theta)
  v = v + sintheta*cross(v,Bhat)+(1-costheta)*cross(cross(v,Bhat),Bhat)
  r = r+v*dt
 #update trail
  trail.append(pos=r)

  KE.plot(pos=(t,1/2*m*mag(v)**2))
