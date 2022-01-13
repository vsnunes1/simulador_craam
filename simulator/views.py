# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.shortcuts import render

from simulator.models import Main, Planet, Star, Orbit, Moon, Spot

# Create your views here.

def home(request):

	return render(request, 'simulator/index.html')

def simulator(request):

	
	
	moons = []

	if int(request.GET['room']) > 0:
		for x in range(1, int(request.GET['room'])+1):

			radius = 'moonRadius_' + str(x)
			mass = 'moonMass_' + str(x)
			albedo = 'moonAlbedo_' + str(x)
			distance = 'moonDistance_' + str(x)
			semiaxis = 'moonSemiaxis_' + str(x)
			period = 'moonPeriod_'+ str(x)
			inclinationAngle = 'moonInclinationAngle_'+ str(x)
			obliquityAngle = 'moonObliquityAngle_'+ str(x)
			eccentricity = 'moonEccentricity_'+ str(x)

			moon = Moon(float(request.GET[radius]), float(request.GET[mass]), float(request.GET[albedo]), float(request.GET[distance]),Orbit(float(request.GET[semiaxis]), float(request.GET[period]), float(request.GET[inclinationAngle]), float(request.GET[obliquityAngle]), float(request.GET[eccentricity])))

			moons.append(moon)


	spots = []

	if int(request.GET['room2']) > 0:
		for x in range(1, int(request.GET['room2'])+1):

			radius = 'spotRadius_' + str(x)
			intensity = 'spotIntensity_' + str(x)
			latitude = 'spotLatitude_' + str(x)
			longitude = 'spotLongitude_' + str(x)
			
			spot = Spot(float(request.GET[radius]), float(request.GET[intensity]), float(request.GET[latitude]), float(request.GET[longitude]))

			spots.append(spot)

	
	star = Star(request.GET['starName'], float(request.GET['starRadius']), float(request.GET['starMass']), float(request.GET['effectiveTemperature']), spots)
	planet = Planet(float(request.GET['planetMass']), float(request.GET['planetRadius']), float(request.GET['atmosphere']), float(request.GET['planetAlbedo']), Orbit(float(request.GET['planetSemiaxis']), float(request.GET['planetPeriod']), float(request.GET['planetInclinationAngle']), float(request.GET['planetObliquityAngle']), float(request.GET['planetEccentricity'])))	
	
	main = Main()
	main.plotImgs(planet, star, moons, request.GET['starColor'])
	main.makeVideo()
	return render(request, 'simulator/simulator.html')

def default(request):
	
	spots = []
	spot = Spot(0.05, 0.5, -30., -55.)
	spots.append(spot)
	spot = Spot(0.06, 0.5, -25., -60.)
	spots.append(spot)
	spot = Spot(0.04, 0.5, -20., -65.)
	spots.append(spot)
	spot = Spot(0.07, 0.5, -15., -70.)
	spots.append(spot)
	spot = Spot(0.1, 0.5, 30., 70.)
	spots.append(spot)
	star = Star("Star", 373., 10, 10, spots)
	planet = Planet(0, 0.066200003, 0, 0, Orbit(19.549999, 9.4341497, 90, 0, 0))
	moon = Moon(0.0288431223213, 0, 0, 4.10784266075, Orbit(0, 0.0751017821823, 90, 0, 0))
	moon2 = Moon(0.0188431223213, 0, 0, 2.10784266075, Orbit(0, 0.0751017821823, 0, 0, 0))

	moons = []
	moons.append(moon)
	moons.append(moon2)	

	main = Main()
	main.plotImgs(planet, star, moons, 'copper')
	main.makeVideo()
	return render(request, 'simulator/simulator.html')

