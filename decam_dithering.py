import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import sys

dit = 20.

cena_ra = (13.+25./60.+27.61507/3600.)*15. 
cena_dec = (-43.-1./60.-8.8053/3600.)

xlim = [cena_ra+7.5/np.cos(cena_dec*np.pi/180.),cena_ra-7.5/np.cos(cena_dec*np.pi/180.)]
ylim = [cena_dec-6.5,cena_dec+5.0]

wcen_ra = (13.+26./60.+47.28/3600.)*15. 
wcen_dec = -47-28./60.-46.1/3600.

gal_coords = [[(13.+27./60.+37.3/3600.)*15.,-41.-28./60.-50./3600.],
			  [(13.+21./60.+47.3/3600.)*15.,-45.-3./60.-42.0/3600.],
			  [(13.+37./60.+39.0/3600.)*15.,-42.-50./60.-49.0/3600.],
			  [(13.+13./60.+11.9/3600.)*15.,-43.-15./60.-56.0/3600.],
			  [(13.+13./60.+9.1/3600.)*15.,-44.-53./60.-24.0/3600.],
			  [(13.+34./60.+47.3/3600.)*15.,-45.-32./60.-51.0/3600.],
			  [(13.+45./60.+0.0/3600.)*15.,-41.-51./60.-40.0/3600.],
			  [(13.+37./60.+25.3/3600.)*15.,-39.-53./60.-48.0/3600.],
			  [(13.+10./60.+32.9/3600.)*15.,-46.-59./60.-27.0/3600.],
			  [(13.+5./60.+2.1/3600.)*15.,-40.-4./60.-58.0/3600.]]

crno_coords = [[(13.+30./60+14.26/3600.)*15.,-41.-53./60.-35.8/3600.],
				[(13.+29./60.+57.34/3600.)*15.,-41.-52./60.-22.6/3600.],
				[(13.+24./60.+10.90/3600.)*15.,-42.-08./60.-24.5/3600.],
			 	[(13.+27./60.+11./3600)*15.,  -42.-29./60.-57./3600.],
			 	[(13.+26./60.+49.16/3600)*15.,-43.-00./60.-01./3600.],
			 	[(13.+27./60.+01.38/3600)*15.,-43.-03./60.-28./3600.],
			 	[(13.+24./60.+10.20/3600)*15.,-43.-23./60.-52./3600.],
			 	[(13.+22./60.+52.7/3600)*15., -43.-32./60.-45./3600.],
			 	[(13.+23./60.+45.4/3600)*15., -43.-39./60.-21.5/3600.],
			 	[(13.+26./60.+28./3600)*15.,  -43.-33./60.-15.6/3600.],
			 	[(13.+27./60.+57.27/3600)*15.,-43.-33./60.-07.8/3600.],
			 	[(13.+28./60.+16.9/3600)*15., -43.-42./60.-18.4/3600.],
			 	[(13.+26./60.+07.42/3600)*15.,-43.-42./60.-28.9/3600.]]
#			 	[13.27:43.23 -41:28:54.6],
#			 	[]
			 	


probes = [[cena_ra-0.1,cena_dec+0.1],[198.0,-45.1],[205.0,-41.0]]
probes = np.array(probes)

#gc_ra, gc_dec = np.loadtxt('wood2007_gc_coords.txt',unpack=True,comments='#')
num, gc_ra, gc_dec = np.loadtxt('CenA_allGC_degcoords.txt',unpack=True,comments='#')

pres = 0.2632

gap_sht = 153. #pixels
gap_lng = 201. #pixels
chip_x = 4096.
chip_y = 2048

gap_sht = gap_sht*pres/3600. #degrees
gap_lng = gap_lng*pres/3600.
chip_x = chip_x*pres/3600.
chip_y = chip_y*pres/3600.

hits_ra = []
hits_dec = []

fileout = open('pointing_coords.txt','w')

def draw_chips(x,y,chip_x0,gap_sht0):
	import matplotlib.patches as patches
	from matplotlib.collections import PatchCollection
	
	row_nums = [2,4,5,6,6,7,7,6,6,5,4,2]
	chips = []
	
	for ii in range(len(row_nums)):
		if ii == 0:
			y0 = (y-6.*chip_y-5.*gap_lng-gap_lng/2.)+11.*chip_y+11.*gap_lng
			chip_x = chip_x0/np.cos(y0*np.pi/180.)
			gap_sht = gap_sht0/np.cos(y0*np.pi/180.)
			x0 = (x-3.*chip_x-3.*gap_sht-chip_x/2.)+4.*chip_x+4.*gap_sht
			for jj in range(row_nums[ii]):
				rect = patches.Rectangle([x0,y0],chip_x,chip_y)
				chips.append(rect)
				hits_ra.append([x0,(x0+chip_x)])
				hits_dec.append([y0,y0+chip_y])
				x0 = x0 - 2.*chip_x - 2.*gap_sht
		if ii == 1:
			y0 = (y-6.*chip_y-5.*gap_lng-gap_lng/2.)+10.*chip_y+10.*gap_lng
			chip_x = chip_x0/np.cos(y0*np.pi/180.)
			gap_sht = gap_sht0/np.cos(y0*np.pi/180.)
			x0 = (x-3.*chip_x-3.*gap_sht-chip_x/2.)+4.*chip_x+5.*gap_sht+chip_x/2.
			for jj in range(row_nums[ii]):
				rect = patches.Rectangle([x0,y0],chip_x,chip_y)
				chips.append(rect)
				hits_ra.append([x0,(x0+chip_x)])
				hits_dec.append([y0,y0+chip_y])
				x0 = x0 - chip_x - gap_sht
		if ii == 2:
			y0 = (y-6.*chip_y-5.*gap_lng-gap_lng/2.)+9.*chip_y+9.*gap_lng
			chip_x = chip_x0/np.cos(y0*np.pi/180.)
			gap_sht = gap_sht0/np.cos(y0*np.pi/180.)
			x0 = (x-3.*chip_x-3.*gap_sht-chip_x/2.)+5.*chip_x+5.*gap_sht
			for jj in range(row_nums[ii]):
				rect = patches.Rectangle([x0,y0],chip_x,chip_y)
				chips.append(rect)
				hits_ra.append([x0,(x0+chip_x)])
				hits_dec.append([y0,y0+chip_y])
				x0 = x0 - chip_x - gap_sht
		if ii == 3:
			y0 = (y-6.*chip_y-5.*gap_lng-gap_lng/2.)+8.*chip_y+8.*gap_lng
			chip_x = chip_x0/np.cos(y0*np.pi/180.)
			gap_sht = gap_sht0/np.cos(y0*np.pi/180.)
			x0 = (x-3.*chip_x-3.*gap_sht-chip_x/2.)+5.*chip_x+6.*gap_sht+chip_x/2.
			for jj in range(row_nums[ii]):
				rect = patches.Rectangle([x0,y0],chip_x,chip_y)
				chips.append(rect)
				hits_ra.append([x0,(x0+chip_x)])
				hits_dec.append([y0,y0+chip_y])
				x0 = x0 - chip_x - gap_sht
		if ii == 4:
			y0 = (y-6.*chip_y-5.*gap_lng-gap_lng/2.)+7.*chip_y+7.*gap_lng
			chip_x = chip_x0/np.cos(y0*np.pi/180.)
			gap_sht = gap_sht0/np.cos(y0*np.pi/180.)
			x0 = (x-3.*chip_x-3.*gap_sht-chip_x/2.)+5.*chip_x+6.*gap_sht+chip_x/2.
			for jj in range(row_nums[ii]):
				rect = patches.Rectangle([x0,y0],chip_x,chip_y)
				chips.append(rect)
				hits_ra.append([x0,(x0+chip_x)])
				hits_dec.append([y0,y0+chip_y])
				x0 = x0 - chip_x - gap_sht
		if ii == 5:
			y0 = (y-6.*chip_y-5.*gap_lng-gap_lng/2.)+6.*chip_y+6.*gap_lng
			chip_x = chip_x0/np.cos(y0*np.pi/180.)
			gap_sht = gap_sht0/np.cos(y0*np.pi/180.)
			x0 = (x-3.*chip_x-3.*gap_sht-chip_x/2.)+6.*chip_x+6.*gap_sht
			for jj in range(row_nums[ii]):
				rect = patches.Rectangle([x0,y0],chip_x,chip_y)
				chips.append(rect)
				hits_ra.append([x0,(x0+chip_x)])
				hits_dec.append([y0,y0+chip_y])
				x0 = x0 - chip_x - gap_sht
		if ii == 6:
			y0 = (y-6.*chip_y-5.*gap_lng-gap_lng/2.)+5.*chip_y+5.*gap_lng
			chip_x = chip_x0/np.cos(y0*np.pi/180.)
			gap_sht = gap_sht0/np.cos(y0*np.pi/180.)
			x0 = (x-3.*chip_x-3.*gap_sht-chip_x/2.)+6.*chip_x+6.*gap_sht
			for jj in range(row_nums[ii]):
				rect = patches.Rectangle([x0,y0],chip_x,chip_y)
				chips.append(rect)
				hits_ra.append([x0,(x0+chip_x)])
				hits_dec.append([y0,y0+chip_y])
				x0 = x0 - chip_x - gap_sht
		if ii == 7:
			y0 = (y-6.*chip_y-5.*gap_lng-gap_lng/2.)+4.*chip_y+4.*gap_lng
			chip_x = chip_x0/np.cos(y0*np.pi/180.)
			gap_sht = gap_sht0/np.cos(y0*np.pi/180.)
			x0 = (x-3.*chip_x-3.*gap_sht-chip_x/2.)+5.*chip_x+6.*gap_sht+chip_x/2.
			for jj in range(row_nums[ii]):
				rect = patches.Rectangle([x0,y0],chip_x,chip_y)
				chips.append(rect)
				hits_ra.append([x0,(x0+chip_x)])
				hits_dec.append([y0,y0+chip_y])
				x0 = x0 - chip_x - gap_sht
		if ii == 8:
			y0 = (y-6.*chip_y-5.*gap_lng-gap_lng/2.)+3.*chip_y+3.*gap_lng
			chip_x = chip_x0/np.cos(y0*np.pi/180.)
			gap_sht = gap_sht0/np.cos(y0*np.pi/180.)
			x0 = (x-3.*chip_x-3.*gap_sht-chip_x/2.)+5.*chip_x+6.*gap_sht+chip_x/2.
			for jj in range(row_nums[ii]):
				rect = patches.Rectangle([x0,y0],chip_x,chip_y)
				chips.append(rect)
				hits_ra.append([x0,(x0+chip_x)])
				hits_dec.append([y0,y0+chip_y])
				x0 = x0 - chip_x - gap_sht
		if ii == 9:
			y0 = (y-6.*chip_y-5.*gap_lng-gap_lng/2.)+2.*chip_y+2.*gap_lng
			chip_x = chip_x0/np.cos(y0*np.pi/180.)
			gap_sht = gap_sht0/np.cos(y0*np.pi/180.)
			x0 = (x-3.*chip_x-3.*gap_sht-chip_x/2.)+5.*chip_x+5.*gap_sht
			for jj in range(row_nums[ii]):
				rect = patches.Rectangle([x0,y0],chip_x,chip_y)
				chips.append(rect)
				hits_ra.append([x0,(x0+chip_x)])
				hits_dec.append([y0,y0+chip_y])
				x0 = x0 - chip_x - gap_sht
		if ii == 10:
			y0 = (y-6.*chip_y-5.*gap_lng-gap_lng/2.)+chip_y+gap_lng
			chip_x = chip_x0/np.cos(y0*np.pi/180.)
			gap_sht = gap_sht0/np.cos(y0*np.pi/180.)
			x0 = (x-3.*chip_x-3.*gap_sht-chip_x/2.)+4.*chip_x+5.*gap_sht+chip_x/2.
			for jj in range(row_nums[ii]):
				rect = patches.Rectangle([x0,y0],chip_x,chip_y)
				chips.append(rect)
				hits_ra.append([x0,(x0+chip_x)])
				hits_dec.append([y0,y0+chip_y])
				x0 = x0 - chip_x - gap_sht
		if ii == 11:
			y0 = (y-6.*chip_y-5.*gap_lng-gap_lng/2.)+0
			chip_x = chip_x0/np.cos(y0*np.pi/180.)
			gap_sht = gap_sht0/np.cos(y0*np.pi/180.)
			x0 = (x-3.*chip_x-3.*gap_sht-chip_x/2.)+4.*chip_x+4.*gap_sht
			for jj in range(row_nums[ii]):
				rect = patches.Rectangle([x0,y0],chip_x,chip_y)
				chips.append(rect)
				hits_ra.append([x0,(x0+chip_x)])
				hits_dec.append([y0,y0+chip_y])
				x0 = x0 - 2.*chip_x - 2.*gap_sht
	return chips

def draw_footprint(point_ra,point_dec,chip_x,gap_sht):
	if sys.argv[1] == 'rect':
		points_x = [point_ra,point_ra+dra,point_ra+dra,point_ra]
		points_y = [point_dec,point_dec,point_dec+ddec,point_dec+ddec]
		p1 = PatchCollection(draw_chips(point_ra,point_dec,chip_x,gap_sht),alpha=0.1)
		p2 = PatchCollection(draw_chips(dra,point_dec,chip_x,gap_sht),alpha=0.1)
		p3 = PatchCollection(draw_chips(dra,ddec,chip_x,gap_sht),alpha=0.1)
		p4 = PatchCollection(draw_chips(point_ra,ddec,chip_x,gap_sht),alpha=0.1)
		collection = [p1,p2,p3,p4]

	if sys.argv[1] == 'rect_centre':
		points_x = [point_ra,point_ra-1.*dra,point_ra+dra,point_ra+dra,point_ra-1.*dra]
		points_y = [point_dec,point_dec-1.*ddec,point_dec-1.*ddec,point_dec+ddec,point_dec+ddec]
		p1 = PatchCollection(draw_chips(point_ra,point_dec,chip_x,gap_sht),alpha=0.05)
		p2 = PatchCollection(draw_chips(point_ra-dra,point_dec-ddec,chip_x,gap_sht),alpha=0.1)
		p3 = PatchCollection(draw_chips(point_ra+dra,point_dec-ddec,chip_x,gap_sht),alpha=0.1)
		p4 = PatchCollection(draw_chips(point_ra+dra,point_dec+ddec,chip_x,gap_sht),alpha=0.1)
		p5 = PatchCollection(draw_chips(point_ra-dra,point_dec+ddec,chip_x,gap_sht),alpha=0.1)
		collection = [p1,p2,p3,p4,p5]

	if sys.argv[1] == 'hexagon':
		points_x = [point_ra,point_ra-1.*dra,point_ra-1.*dra,point_ra,dra,dra]
		points_y = [point_dec,point_dec+ddec,point_dec+2.*ddec,point_dec+3.*ddec,point_dec+2.*ddec,point_dec+ddec]
		p1 = PatchCollection(draw_chips(point_ra,point_dec,chip_x,gap_sht),alpha=0.1)
		p2 = PatchCollection(draw_chips(point_ra-dra,point_dec+ddec,chip_x,gap_sht),alpha=0.1)
		p3 = PatchCollection(draw_chips(point_ra-dra,point_dec+2.*ddec,chip_x,gap_sht),alpha=0.1)
		p4 = PatchCollection(draw_chips(point_ra,point_dec+3.*ddec,chip_x,gap_sht),alpha=0.1)
		p5 = PatchCollection(draw_chips(point_ra+dra,point_dec+2.*ddec,chip_x,gap_sht),alpha=0.1)
		p6 = PatchCollection(draw_chips(point_ra+dra,point_dec+ddec,chip_x,gap_sht),alpha=0.1)
		collection = [p1,p2,p3,p4,p5,p6]
	
	return points_x, points_y, collection

def get_pointings(cena_ra,cena_dec,chip_x,gap_sht,chip_y,gap_lng,type):
	if type == 'orig':
		pointings = [[cena_ra,cena_dec],
[cena_ra+(3.*chip_x+3.*gap_sht+chip_x/2.)/np.cos((cena_dec+6.*chip_y+6.*gap_lng)*np.pi/180.),cena_dec+6.*chip_y+6.*gap_lng],
			 [cena_ra,cena_dec+10.*chip_y+10.*gap_lng],
			 [cena_ra-(3.*chip_x+3.*gap_sht+chip_x/2.)/np.cos((cena_dec+6.*chip_y+6.*gap_lng)*np.pi/180.),cena_dec+6.*chip_y+6.*gap_lng],
			 [cena_ra-(7.*chip_x+7.*gap_sht)/np.cos(cena_dec*np.pi/180.),cena_dec],
			 [cena_ra-(3.*chip_x+3.*gap_sht+chip_x/2.)/np.cos((cena_dec-6.*chip_y-6.*gap_lng)*np.pi/180.),cena_dec-6.*chip_y-6.*gap_lng],
			 [cena_ra,cena_dec-10.*chip_y-10.*gap_lng],
			 [cena_ra+(3.*chip_x+3.*gap_sht+chip_x/2.)/np.cos((cena_dec-6.*chip_y-6.*gap_lng)*np.pi/180.),cena_dec-6.*chip_y-6.*gap_lng],
			 [cena_ra+(7.*chip_x+7.*gap_sht)/np.cos(cena_dec*np.pi/180.),cena_dec],
			 #Outside bottom circle arc
			 [cena_ra+(10.*chip_x+10.*gap_sht)/np.cos((cena_dec-6.*chip_y-6.*gap_lng)*np.pi/180.),cena_dec-6.*chip_y-6.*gap_lng],
			 [cena_ra+(7.*chip_x+7.*gap_sht)/np.cos((cena_dec-10.*chip_y-10.*gap_lng)*np.pi/180.),cena_dec-10.*chip_y-10.*gap_lng],
			 [cena_ra+(3.*chip_x+3.*gap_sht+chip_x/2.)/np.cos((cena_dec-16.*chip_y-16.*gap_lng)*np.pi/180.),cena_dec-16.*chip_y-16.*gap_lng],
			 [cena_ra,cena_dec-20.*chip_y-20.*gap_lng],
			 [cena_ra-(3.*chip_x+3.*gap_sht+chip_x/2.)/np.cos((cena_dec-16.*chip_y-16.*gap_lng)*np.pi/180.),cena_dec-16.*chip_y-16.*gap_lng],
			 [cena_ra-(7.*chip_x+7.*gap_sht)/np.cos((cena_dec-10.*chip_y-10.*gap_lng)*np.pi/180.),cena_dec-10.*chip_y-10.*gap_lng],
			 [cena_ra-(10.*chip_x+10.*gap_sht)/np.cos((cena_dec-6.*chip_y-6.*gap_lng)*np.pi/180.),cena_dec-6.*chip_y-6.*gap_lng],
 			 #Outside top circle arc
			 [cena_ra-(10.*chip_x+10.*gap_sht)/np.cos((cena_dec+6.*chip_y+6.*gap_lng)*np.pi/180.),cena_dec+6.*chip_y+6.*gap_lng],
			 [cena_ra-(7.*chip_x+7.*gap_sht)/np.cos((cena_dec+10.*chip_y+10.*gap_lng)*np.pi/180.),cena_dec+10.*chip_y+10.*gap_lng],
			 [cena_ra-(3.*chip_x+3.*gap_sht+chip_x/2.)/np.cos((cena_dec+16.*chip_y+16.*gap_lng)*np.pi/180.),cena_dec+16.*chip_y+16.*gap_lng],
			[cena_ra,cena_dec+20.*chip_y+20.*gap_lng],
			 [cena_ra+(3.*chip_x+3.*gap_sht+chip_x/2.)/np.cos((cena_dec+16.*chip_y+16.*gap_lng)*np.pi/180.),cena_dec+16.*chip_y+16.*gap_lng],
			[cena_ra+(7.*chip_x+7.*gap_sht)/np.cos((cena_dec+10.*chip_y+10.*gap_lng)*np.pi/180.),cena_dec+10.*chip_y+10.*gap_lng],
			[cena_ra+(10.*chip_x+10.*gap_sht)/np.cos((cena_dec+6.*chip_y+6.*gap_lng)*np.pi/180.),cena_dec+6.*chip_y+6.*gap_lng]
			]
	if type == 'test':
		pointings = [[cena_ra,cena_dec],
			 [cena_ra+(4.*chip_x-0.5*gap_sht+chip_x/2.)/np.cos((cena_dec+7.*chip_y+7.*gap_lng)*np.pi/180.),cena_dec+7.*chip_y+7.*gap_lng],
			 [cena_ra-(chip_x+1.5*gap_sht)/np.cos((cena_dec+11.*chip_y+11.*gap_lng)*np.pi/180.),cena_dec+11.*chip_y+11.*gap_lng],
			 [cena_ra-(5.*chip_x+0.5*gap_sht+chip_x/2.)/np.cos((cena_dec+4.*chip_y+4.*gap_lng)*np.pi/180.),cena_dec+4.*chip_y+4.*gap_lng],
			 [cena_ra-(4.*chip_x-4.5*gap_sht+chip_x/2.)/np.cos((cena_dec-7.*chip_y-7.*gap_lng)*np.pi/180.),cena_dec-7.*chip_y-7.*gap_lng],
			 [cena_ra+(chip_x+gap_sht)/np.cos((cena_dec-11.*chip_y-11.*gap_lng)*np.pi/180.),cena_dec-11.*chip_y-11.*gap_lng],
			 [cena_ra+(5.*chip_x-3.5*gap_sht+chip_x/2.)/np.cos((cena_dec-4.*chip_y-4.*gap_lng)*np.pi/180.),cena_dec-4.*chip_y-4.*gap_lng],
			 [cena_ra+(10.*chip_x-3.0*gap_sht)/np.cos((cena_dec+3.*chip_y+3.*gap_lng)*np.pi/180.),cena_dec+3.*chip_y+3.*gap_lng],
			 [cena_ra+(10.*chip_x+3.*gap_sht+chip_x/2.)/np.cos((cena_dec-8.*chip_y-8.*gap_lng)*np.pi/180.),cena_dec-8.*chip_y-8.*gap_lng],
			 [cena_ra+(6.*chip_x-7.5*gap_sht+chip_x/2.)/np.cos((cena_dec-15.*chip_y-15.*gap_lng)*np.pi/180.),cena_dec-15.*chip_y-15.*gap_lng],
			 [cena_ra+(2.*chip_x+gap_sht)/np.cos((cena_dec-22.*chip_y-22.*gap_lng)*np.pi/180.),cena_dec-22.*chip_y-22.*gap_lng],
			 [cena_ra-(3.*chip_x-8.5*gap_sht+chip_x/2.)/np.cos((cena_dec-18.*chip_y-18.*gap_lng)*np.pi/180.),cena_dec-18.*chip_y-18.*gap_lng],
			 [cena_ra-(8.*chip_x-21.*gap_sht)/np.cos((cena_dec-24.*chip_y-24.*gap_lng)*np.pi/180.),cena_dec-24.*chip_y-24.*gap_lng],
			 [cena_ra-(9.*chip_x-13.*gap_sht)/np.cos((cena_dec-13.*chip_y-13.*gap_lng)*np.pi/180.),cena_dec-13.*chip_y-13.*gap_lng],
			 [cena_ra-(10.*chip_x-4.5*gap_sht)/np.cos((cena_dec-2.*chip_y-2.*gap_lng)*np.pi/180.),cena_dec-2.*chip_y-2.*gap_lng],
			 [cena_ra-(11.*chip_x+4.5*gap_sht)/np.cos((cena_dec+9.*chip_y+9.*gap_lng)*np.pi/180.),cena_dec+9.*chip_y+9.*gap_lng],
			 [cena_ra-(12.*chip_x+14.*gap_sht)/np.cos((cena_dec+20.*chip_y+20.*gap_lng)*np.pi/180.),cena_dec+20.*chip_y+20.*gap_lng],
			 [cena_ra-(6.*chip_x+6.*gap_sht+chip_x/2.)/np.cos((cena_dec+15.*chip_y+15.*gap_lng)*np.pi/180.),cena_dec+15.*chip_y+15.*gap_lng],
			 [cena_ra-(2.*chip_x+4.*gap_sht)/np.cos((cena_dec+22.*chip_y+22.*gap_lng)*np.pi/180.),cena_dec+22.*chip_y+22.*gap_lng],
			 [cena_ra+(3.*chip_x+1.*gap_sht+chip_x/2.)/np.cos((cena_dec+18.*chip_y+18.*gap_lng)*np.pi/180.),cena_dec+18.*chip_y+18.*gap_lng],
			 [cena_ra+(9.*chip_x+3.*gap_sht)/np.cos((cena_dec+14.*chip_y+14.*gap_lng)*np.pi/180.),cena_dec+14.*chip_y+14.*gap_lng],
			 [cena_ra+(14.*chip_x+3.5*gap_sht+chip_x/2.)/np.cos((cena_dec+9.*chip_y+9.*gap_lng)*np.pi/180.),cena_dec+9.*chip_y+9.*gap_lng],
			 [cena_ra+(8.*chip_x+8.*gap_sht)/np.cos((cena_dec+25.*chip_y+25.*gap_lng)*np.pi/180.),cena_dec+25.*chip_y+25.*gap_lng],
			 [cena_ra-(0.5*chip_x-12.0*gap_sht)/np.cos((cena_dec-29.*chip_y-29.*gap_lng)*np.pi/180.),cena_dec-29.*chip_y-29.*gap_lng],
#			 [cena_ra+(3.*chip_x-0.5*gap_sht)/np.cos((cena_dec-33.*chip_y-33.*gap_lng)*np.pi/180.),cena_dec-33.*chip_y-33.*gap_lng],
			]
	return pointings


fig, ax = plt.subplots(figsize=(12,10))

#dra = 60./3600.
#ddec = 75./3600.

dra = 41./3600.
ddec = 54.0/3600.

#print dra*3600., ddec*3600.
#exit()

pointings = get_pointings(cena_ra,cena_dec,chip_x,gap_sht,chip_y,gap_lng,'test')
pointings_all_dec = []
for ii in range(len(pointings)):
	points = [[pointings[ii][0],pointings[ii][1]],
			  [pointings[ii][0]-1.*dra,pointings[ii][1]-1.*ddec],
			  [pointings[ii][0]+dra,pointings[ii][1]-1.*ddec],
			  [pointings[ii][0]+dra,pointings[ii][1]+ddec],
			  [pointings[ii][0]-1.*dra,pointings[ii][1]+ddec]]
	print >> fileout, "p%i" % (ii+1)
	for jj in range(len(points)):
		print >> fileout, points[jj][0], points[jj][1]

fileout.close()

for ii in range(len(pointings)):
	point_ra = pointings[ii][0]
	point_dec = pointings[ii][1]	

	print point_ra, point_dec

	points_x, points_y, collection = draw_footprint(point_ra,point_dec,chip_x,gap_sht)

	for patches in collection:
		ax.add_collection(patches)

	ax.plot(points_x,points_y,'ko')
	ax.annotate('%i' % (ii+1),(pointings[ii][0]-0.1,pointings[ii][1]+0.1))

print "Total hits = %i" % len(hits_ra)
hits_ra1 = []
hits_ra0 = []
hits_dec0 = []
hits_dec1 = []
for ii in range(len(hits_ra)):
	hits_ra0.append(hits_ra[ii][0])
	hits_ra1.append(hits_ra[ii][1])
	hits_dec0.append(hits_dec[ii][0])
	hits_dec1.append(hits_dec[ii][1])

pcheck_ra = np.linspace(np.min(hits_ra0),np.max(hits_ra1),1e4)
pcheck_dec = np.linspace(np.min(hits_dec0),np.max(hits_dec1),1e4)

counts = []

if sys.argv[2] == 'probe':
	print "Probing coverage..."
	for ii in range(len(probes)):
		print "Probe %i, RA...." % (ii+1)
		cnt_ra = []	
		for ra in pcheck_ra:
			mask = np.where(hits_dec0 <= probes[ii][1],1,0)
			chk_ra0 = np.compress(mask,hits_ra0)
			chk_ra1 = np.compress(mask,hits_ra1)
			chk_dec0 = np.compress(mask,hits_dec0)
			chk_dec1 = np.compress(mask,hits_dec1)
			mask = np.where(chk_dec1 >= probes[ii][1],1,0)
			chk_ra0 = np.compress(mask,chk_ra0)
			chk_ra1 = np.compress(mask,chk_ra1)
			chk_dec0 = np.compress(mask,chk_dec0)
			chk_dec1 = np.compress(mask,chk_dec1)
			mask = np.where(chk_ra0 <= ra,1,0)
			chk_ra0 = np.compress(mask,chk_ra0)
			chk_ra1 = np.compress(mask,chk_ra1)
			chk_dec0 = np.compress(mask,chk_dec0)
			chk_dec1 = np.compress(mask,chk_dec1)
			mask = np.where(chk_ra1 >= ra,1,0)
			chk_ra0 = np.compress(mask,chk_ra0)
			chk_ra1 = np.compress(mask,chk_ra1)
			chk_dec0 = np.compress(mask,chk_dec0)
			chk_dec1 = np.compress(mask,chk_dec1)
		
			cnt_ra.append(len(chk_ra0))
	
		print "Probe %i, Dec...." % (ii+1)
		cnt_dec = []	
		for dec in pcheck_dec:
			mask = np.where(hits_ra0 <= probes[ii][0],1,0)
			chk_ra0 = np.compress(mask,hits_ra0)
			chk_ra1 = np.compress(mask,hits_ra1)
			chk_dec0 = np.compress(mask,hits_dec0)
			chk_dec1 = np.compress(mask,hits_dec1)
			mask = np.where(chk_ra1 >= probes[ii][0],1,0)
			chk_ra0 = np.compress(mask,chk_ra0)
			chk_ra1 = np.compress(mask,chk_ra1)
			chk_dec0 = np.compress(mask,chk_dec0)
			chk_dec1 = np.compress(mask,chk_dec1)
			mask = np.where(chk_dec0 <= dec,1,0)
			chk_ra0 = np.compress(mask,chk_ra0)
			chk_ra1 = np.compress(mask,chk_ra1)
			chk_dec0 = np.compress(mask,chk_dec0)
			chk_dec1 = np.compress(mask,chk_dec1)
			mask = np.where(chk_dec1 >= dec,1,0)
			chk_ra0 = np.compress(mask,chk_ra0)
			chk_ra1 = np.compress(mask,chk_ra1)
			chk_dec0 = np.compress(mask,chk_dec0)
			chk_dec1 = np.compress(mask,chk_dec1)
		
			cnt_dec.append(len(chk_ra0))
	
		counts.append([cnt_ra,cnt_dec])
	
	counts = np.array(counts)

theta = np.linspace(0,2*np.pi,1e5)
r = np.arctan(3.e5/3.8e6)*180./np.pi
y = cena_dec+r*np.sin(theta)
x = cena_ra-r*np.cos(theta)/np.cos(y*np.pi/180.)
print r, cena_ra, cena_dec
print np.min(x), np.max(x)
print np.min(y), np.max(y)

ax.plot(x,y,'r--')

wcen, = ax.plot([wcen_ra],[wcen_dec],'o',color='blue',ms=20)

#for ii in range(len(probes)):
#	ax.plot([xlim[0],xlim[1]],[probes[ii][1],probes[ii][1]],'k--')
#	ax.plot([probes[ii][0],probes[ii][0]],[ylim[0],ylim[1]],'k--')
#	ax.plot([probes[ii][0]],[probes[ii][1]],'kx',ms=20)

for ii in range(len(gal_coords)):
	dwarf, = ax.plot([gal_coords[ii][0]],[gal_coords[ii][1]],'^',color='darkorange',ms=15)

for ii in range(len(crno_coords)):
	crno, = ax.plot([crno_coords[ii][0]],[crno_coords[ii][1]],'o',color='maroon',ms=12)

gcs, = ax.plot(gc_ra,gc_dec,'g,')

cena, = ax.plot([cena_ra],[cena_dec],'*',color='green',ms=25)

ax.legend((cena,gcs,wcen,dwarf,crno),(r'NGC$\,$5128',r'Confirmed GCs $$',r'$\omega$ Centauri',r'Known Dwarfs $$',r'Dwarf Candidates $$'),loc=3,numpoints=1)

ax.set_xlim(xlim[0],xlim[1])
ax.set_ylim(ylim[0],ylim[1])
# ax.set_xlim(208.8,195.0)
# ax.set_ylim(ylim[0],ylim[1])
ax.tick_params(axis='both',labelsize=18)
ax.set_xlabel(r'$\alpha$ [deg]',fontsize=20)
ax.set_ylabel(r'$\delta$ [deg]',fontsize=20)

plt.savefig('dithering_pattern.pdf')
#plt.savefig('dithering_pattern_4poster.pdf')

if sys.argv[2] == 'probe':
	plt.figure(figsize=(10,10))
	cnt = 0
	for ii in range(len(probes)):
		cnt = cnt + 1
		plt.subplot('%i%i%i' % (len(probes),2,cnt))
		plt.plot(pcheck_ra,counts[ii][0]*dit,'b')
		plt.xlim(xlim[0],xlim[1])
		plt.xlabel(r'RA [deg]')
		plt.ylabel(r'Exp. Time [s]')
	
		cnt = cnt + 1
		plt.subplot('%i%i%i' % (len(probes),2,cnt))
		plt.plot(pcheck_dec,counts[ii][1]*dit,'b')
		plt.xlim(ylim[0],ylim[1])
		plt.xlabel(r'Dec [deg]')
		plt.ylabel(r'Exp. Time [s]')
	
	plt.savefig('decam_probes.pdf')

plt.show()
	