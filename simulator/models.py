from __future__ import unicode_literals
from django.db import models
from SMOOTH import blur_image, smooth
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math
import numpy as np
import os
import sys


class Star (object):
    
    #MATRIX LEN
    matrixLen = 856
    
    #Maximus middle of star intensity
    maximusIntensity=240.

    #Input
    coeficienteHum = 0.5
    coeficienteDois = 0.3

    imgStar = np.asarray([])
    
    def __init__(self, name, radius, mass, effectiveTemperature, spots):

        self.name = name
        self.radius = radius
        self.mass = mass
        self.effectiveTemperature = effectiveTemperature
        self.spots = spots
        self.ImgStar()
        self.createSpots()
        
        
    def ImgStar(self):
    	global imgStar
    
        star = [[ 0.0 for i in range(self.matrixLen)] for j in range(self.matrixLen)]
        
        for j in range(len(star)):
            for i in range(len(star[j])):
                
                #Distance to de Middle of Star
                middleDistance = math.sqrt(pow(i-self.matrixLen/2,2) + pow(j-self.matrixLen/2,2))
                
                if middleDistance <= self.radius:
                    cosTheta = math.sqrt(1-pow(middleDistance/self.radius,2))
                    star[i][j] = int(self.maximusIntensity * (1 - self.coeficienteHum *(1 - cosTheta) - self.coeficienteDois * (pow(1 - cosTheta, 2))))
        
        imgStar = np.asarray(star)
        return np.asarray(star)

    def returnimgStar(self):
    	global imgStar
    	return imgStar

    def createSpots(self):
		global imgStar


		for i in range(0, len(self.spots)):

			spotTeste = self.spots[i]

			# Parametros da mancha
			raioMancha = self.radius*0.05  # raio em funcao do raio da estrela em pixels

			estrela = self.returnimgStar()


			#######  inserir manchas

			#posicao da mancha em pixels em relacao ao centro da estrela
			ys=self.radius*np.sin(spotTeste.latitude)  
			xs=self.radius*np.cos(spotTeste.latitude)*np.sin(spotTeste.longitude)
			anguloHelio=np.arccos(np.cos(spotTeste.latitude)*np.cos(spotTeste.longitude))

			Ny, Nx = self.returnimgStar().shape

			# efeito de projecao pela mancha estar a um anguloHeliocentrico do centro da estrela - elipcidade
			yy = ys + Ny/2 # posicao em pixel com relacao a origem da matriz
			xx = xs + Nx/2 # posicao em pixel com relacao a origem da matriz

			kk = np.arange(Ny * Nx)
			vx = kk-Nx*np.int64(1.*kk/Nx) - xx
			vy = kk/Ny - yy

			# angulo de rotacao da mancha
			anguloRot=np.abs(np.arctan(ys/xs))    # em radianos
			if spotTeste.latitude*spotTeste.longitude > 0: anguloRot=-anguloRot

			ii, = np.where((((vx*np.cos(anguloRot)-vy*np.sin(anguloRot))/np.cos(anguloHelio))**2+(vx*np.sin(anguloRot)+vy*np.cos(anguloRot))**2) < raioMancha**2)
			 
			spot = np.zeros(Ny * Nx) + 1
			        
			spot[ii]=spotTeste.intensity
			spot = spot.reshape([Ny, Nx])

			imgStar = estrela * spot


class Planet (object):
    
    def __init__(self, mass, radius, atmosphere, albedo, orbit):
        self.mass = mass
        self.radius = radius #(in Rstar)
        self.atmosphere = atmosphere
        self.albedo = albedo
        self.orbit = orbit
        
    def tlstep(self, Rs_pix, ImgStar):
        Rpl_pix = self.radius * Rs_pix
        Rorb_pix = self.orbit.semiaxis * Rs_pix
        
        tdur = self.orbit.period*24/(2*np.pi * Rorb_pix) * (Rs_pix + Rpl_pix) # transit duration in hours / 2pi
        
        tdur = self.orbit.period*24/(2*np.pi * Rorb_pix) * (Rs_pix + Rpl_pix) # transit duration in hours / 2pi

        nstep = int(tdur*60.*4/self.orbit.returndt())# + 0.5) # number of steps
        if nstep%2!=0: nstep = nstep+1
        tstep = np.arange(-nstep/2., nstep/2., 1) * self.orbit.returndt()/60. # in hours
        Lstep = np.zeros(nstep) + ImgStar.sum()# / ImgStar.sum()
        
        tlstep = []
        tlstep.append(tstep)
        tlstep.append(Lstep)
        
        return tlstep
    
    def planetOrbit(self,Rs_pix, ImgStar):
        tstep = self.tlstep(Rs_pix, ImgStar)
        theta_p = 2*np.pi * tstep[0] / (self.orbit.period*24.) - np.pi/2.
        xp = self.orbit.semiaxis * Rs_pix * np.cos(theta_p)
        yp = self.orbit.semiaxis * Rs_pix * np.sin(theta_p) * np.cos(self.orbit.inclinationAngle)
        
        pair = []
        pair.append(xp)
        pair.append(yp)
        return pair
    
    def createPlanet(self, Ny, Nx, kk, y0, x0, Rpl_pix, Rs_pix, moon):
        planet = np.zeros(Ny * Nx) + 1
        # Planet:
        ii, = np.where((kk/Nx-y0)**2+(kk-Nx*np.int64(1.*kk/Nx)-x0)**2 < Rpl_pix**2)
        planet[ii] = 0.
        
        if len(moon) > 0:            
            for x in Range(0, len(moon)):
                moon[x]

class Moon (object):
    
    #pm = 0.0751017821823 #moon period
    #rm = 0.0288431223213 # moon radius
    #dm = 4.10784266075 # moon distance
    tm0 = 1.15612181491 # moon first transit time
    #kind = 'big-fst_'
    
    pos = np.random.choice([-1, 1])
    
    
    def __init__(self, radius, mass, albedo, distance, orbit):
        self.radius = radius
        self.mass = mass
        self.albedo = albedo
        self.distance = distance
        self.orbit = orbit
        
        
        
    # moon orbit in equatorial plane of planet
    def moonOrbit(self, Rs_pix, Rpl_pix, tstep):
        self.distance = self.distance * self.pos    
        Rmoon = self.radius * Rs_pix
        dmoon = self.distance * Rpl_pix
        theta_m0 = self.tm0
        
        theta_m = 2*np.pi * tstep / (self.orbit.period*24.) - theta_m0
        xm = dmoon * np.cos(theta_m)
        ym = dmoon * np.sin(theta_m) * np.cos(self.orbit.inclinationAngle) 
        
        pair = []
        pair.append(xm)
        pair.append(ym)
        
        return pair
    
    def dMoon(self, Rpl_pix):
        return self.distance * Rpl_pix
    
    def rMoon(self, Rs_pix):
        return self.radius * Rs_pix
        
class Orbit (object):
    def __init__(self, semiaxis, period, inclinationAngle, obliquityAngle, eccentricity):
        self.semiaxis = semiaxis #(in Rstar)
        self.period = period #(in days), assumed circular orbit
        self.inclinationAngle = np.radians(inclinationAngle) #(in rad)
        self.obliquityAngle = obliquityAngle
        self.eccentricity = eccentricity
        
        
    def returndt(self): #intervalo entre os pontos
        return 2 #ou 20 testar

class Spot (object):

    def __init__(self, radius, intensity, latitude, longitude):
        self.radius = radius
        self.intensity = intensity
        self.latitude = np.radians(latitude)
        self.longitude = np.radians(longitude)

class Main (object):

	def __init__(self):
		print 'Initialize'

	def makeVideo(self):
		try:
			os.remove("simulator/static/video/test.mp4")
		except OSError:
			pass
		
		os.system("avconv -r 10 -start_number 40 -i simulator/output/frames/img_%d.png -b:v 1000k simulator/static/video/test.mp4")

		
	def plotImgs(self, planet, star, moons, maps):
	
		plt.rcParams['xtick.minor.visible'] = True
		plt.rcParams['ytick.minor.visible'] = True
		plt.rcParams['legend.borderaxespad'] = .3
		plt.rcParams['legend.handlelength'] = 1.
		
		sMa = planet.orbit.semiaxis #semi-major axis of planetary orbit (in Rstar)
		Rs_pix = star.radius #stellar radius (in ~pixel?)
		Rplan = planet.radius #planetary radius (in Rstar) (1 Rjup = 6.9911e4 km = 0.100447 Rsun)
		ImgStar = star.returnimgStar()

		Ny, Nx = ImgStar.shape

		Rpl_pix = Rplan * Rs_pix

		coordp = planet.planetOrbit(Rs_pix, ImgStar)
		xp = coordp[0]
		yp = coordp[1]

		tlstep = planet.tlstep(Rs_pix, ImgStar)
		tstep = tlstep[0]
		Lstep = tlstep[1]
		
		
		#VARIAVEIS DE MATRIZ
		kk = np.arange(Ny * Nx) 

		q, = np.where((xp+(Nx/2) + 98 > 0) & (xp+(Nx/2) - 98 < Nx))
		i = int(np.mean(q))
		x0 = xp[i] + Nx/2
		y0 = yp[i] + Ny/2
		planet = np.zeros(Ny * Nx) + 1
		
		# Planet:
		ii, = np.where((kk/Nx-y0)**2+(kk-Nx*np.int64(1.*kk/Nx)-x0)**2 < Rpl_pix**2)
		planet[ii] = 0.

		
		# Moons:
		if len(moons) > 0:
			for x in range(0, len(moons)):
				
				ll, = np.where((kk/Nx-(y0-moons[x].moonOrbit(Rs_pix, Rpl_pix, tstep)[1][i]))**2+(kk-Nx*np.int64(kk/Nx)-(x0-moons[x].moonOrbit(Rs_pix, Rpl_pix, tstep)[0][i]))**2 < moons[x].rMoon(Rs_pix)**2)
				planet[ll]=0.
		
		planet = planet.reshape([Ny, Nx])
	            
		Lmin = np.sum(ImgStar*planet) / ImgStar.sum()

		
		color = '.9'
		g = 1

		#print 'Loop start. '
		for i in q:
			planet = np.zeros(Ny * Nx) + 1

			x0 = xp[i] + Nx/2
			y0 = yp[i] + Ny/2

			# Planet:
			ii, = np.where((kk/Nx-y0)**2+(kk-Nx*np.int64(1.*kk/Nx)-x0)**2 < Rpl_pix**2)
			planet[ii] = 0.

			# Moon:
			if len(moons) > 0:
				for x in range(0, len(moons)):
					ll, = np.where((kk/Nx-(y0-moons[x].moonOrbit(Rs_pix, Rpl_pix, tstep)[1][i]))**2+(kk-Nx*np.int64(kk/Nx)-(x0-moons[x].moonOrbit(Rs_pix, Rpl_pix, tstep)[0][i]))**2 < moons[x].rMoon(Rs_pix)**2)
					planet[ll]=0.
			planet = planet.reshape([Ny, Nx])

			# Smooth image
			img = blur_image(ImgStar*planet, 5)

			Lstep[i] = np.sum(img)
	        #if L!=0: Lstep[i] = np.sum(img)* np.random.np.random.normal(1., 3e-4)
	        #else:    Lstep[i] = np.sum(img)
			lc = Lstep / ImgStar.sum()
	        #lc = lc + (1 - np.nanmax(lc))

	        # Plotting ============================================================
			fig = plt.figure(figsize=(8, 6), dpi=128)#, facecolor='k', edgecolor='k')
			gs = gridspec.GridSpec(4, 1)#, width_ratios=[1, 1])
			#ax = fig.add_subplot(111)
			ax = plt.subplot(gs[:3, :])
			ax.imshow(img, cmap=maps, interpolation='gaussian', origin='lower')
			ax = plt.subplot(gs[3:, :])
			ax.plot(tstep[:i], lc[:i], '.-', lw=.8, ms=3, c='#ef9f4f')
			ax.set_xlim(tstep[0], tstep[-1])#*1.1)
			ax.set_ylim(Lmin*.999, 1.001)
			ax.set_xlabel('time (h)', fontweight='bold')
			ax.set_ylabel(r'flux (L$\star$)', fontweight='bold')
			ax.grid(which='major', c=color, alpha=.6, lw=.6)
			ax.grid(which='minor', c=color, alpha=.3, lw=.3)
			ax.set_facecolor('k')
			for sp in ('left','bottom', 'right', 'top'):
				ax.spines[sp].set_color(color)
			ax.xaxis.label.set_color(color)
			ax.yaxis.label.set_color(color)
			ax.tick_params(axis='both', which='both', colors=color)
			fig.tight_layout()
			fig.savefig('simulator/output/frames/img_' + str(i) + '.png', dpi=128, facecolor='k')
			plt.close('all')
			g += 1
	        # =====================================================================
		print 'Process Done!'
		return "Ok!"


