# -*- coding: mbcs -*-

import os
import os.path
import time
import math
import numpy as np
from abaqus import *
from abaqusConstants import *
from caeModules import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from driverUtils import *
import visualization
import odbAccess
import regionToolset
import csv
import xyPlot
import displayGroupOdbToolset as dgo
import __main__
import time
import shutil
import subprocess
csv.list_dialects()

executeOnCaeStartup()

################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################
''' Geometrically and materially nonlinear analysis of timber beam-columns '''
################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################

# A python script for automized and parameterised analysis of timber beam-columns with Abaqus.

# Structure of the script
# 	84	Input		122		Geometry
#					333		Boundary Conditions & Loading	
#					399		Material						
#					456		Analysis						
#					477		Saving the model
#					485		Material data
# 	1004 Computation of further input parameters -> from here on the code does not have to be changed any more
#				    1005	Beam heights
# 	1026 Modelling 	1028	Parts
#					1122	Sections
#					1247	Material Orientation
#					1282	Assembly
#					1292	Reference nodes
#					1563	Sets
#					2248	Boundary Conditions
#					2337	Fricton at supports 
#					2367	Springs of the bracing
#					2429	Coupling
#					2553	Loading
#					2818	Output settings
#	2855	Create & Submit Job
# 	2879	Postprocessing
#					2881	Output
#					2898	Fin
#
# Tips:
# - define numbers as integers, e.g. L = 8000. instead of L = 8000

error_exit = 0							# Pointer of fatal error_exitoccurs and analysis has to be aborted
################################################################################################################################################################################################################################################
''' Input '''
################################################################################################################################################################################################################################################
safe_model = 1						# Safe the model, 0 = no, 1 = yes
submit_job = 1						# Submit the job, 0 = no, 1 = yes
del_job = 0							# Delete job after computation, 0 = no, 1 = yes
del_files = 0						# Delete files in work_directory after computation, 0 = no, 1 = yes

# Variation of input parameters
type_var = 'manu'					# Manual = 'manu', automatic = 'auto', when using a KaReMo 'auto' needs to be picked

# Manual variation of input
# L_var = [2000.,2000.,2000.,3000.,3000.,3000.,4000.,4000.,4000.]			# Variation needed when activating loop in line 87

# Automatic variation of input - here material properties and imperfections
nS = 200							# Size of sample for Monte Carlo Simulation
scatter_geomImp_y = 1				# Consider scattering geometrical imperfections according to scatter_geomImp.txt in y direction, 0 = no, 1 = yes
scatter_geomImp_z = 0				# Consider scattering geometrical imperfections according to scatter_geomImp.txt in z direction, 0 = no, 1 = yes
KaReMo = 0							# Consider scattering material values according to KaReMo, 0 = no, 1 = yes; other material values ar picked according to the chosen material and the material data table
KaReMo_v = 0						# Consider scattering material values according to KaReMo_v, 0 = no, 1 = yes; other material values ar picked according to the chosen material and the material data table
KaReMo_y = 0						# Consider scattering material values according to KaReMo_y, 0 = no, 1 = yes; other material values ar picked according to the chosen material and the material data table
Hlam = 40							# Height of lamellae for scattering material values, is also used as element height
newLam = 1							# Start each lamellae with new board
noLTB = 0							# Prevent LTB with horizontal supports at all nodes, 0 = no, 1 = yes
noIPB = 0							# Prevent in-plane buckling with horizontal supports at all nodes, 0 = no, 1 = yes
# Read scattering imperfections
if type_var == 'auto' and (scatter_geomImp_y == 1 or scatter_geomImp_z == 1):
	data_Imp = np.zeros((nS,4))
	c = 1
	i = 0
	with open('data_Imp.txt') as f:
		for x in f:
			if c >= 2 and c < 2+nS:
				x = x.split(',')
				data_Imp[i] = ([float(x[0]),float(x[1]),float(x[2]),float(x[3])])
				i = i + 1
			c = c + 1

# Loop for variation of input parameters
for i in range(1): 				# Computation of one beam
# for i in range(len(H_var)): 		# Computation of multiple beams with varying properties
# for i in range(nS):				# Monte Carlo Simulation with scattering geometric imperfections and/or material properties

	########################################
	''' Geometry '''
	########################################
	
	''' Cross-section Dimensions '''
	L 				= 3450.			# Member lenght [mm], alternatively: L = L_var[i]
	H				= 4*120.		# Member height [mm]
	B				= 120.			# Member width [mm]
	
	''' Meshing '''					# Always chose even numbers
	nEL			    = int(L/150)	# Number of elements over the length / x direction, 2 <= nEL <= 98
	nEB				= 6				# Number of elements over the width  / z direction, 2 <= nEB <= 98
	nEH				= int(H/Hlam)	# Number of elements over the height / y direction, 4 <= nEH <= 98
	sizeEL 			= L/float(nEL)	# Element size in length / x direction	
	sizeEB 			= B/float(nEB)	# Element size in width  / y direction
	sizeEH 			= H/float(nEH)	# Element size in height / z direction	
	nNL	 			= nEL + 1		# Number of nodes in length / x direction
	nNH				= nEH + 1		# Number of nodes in height / z direction
	nNB				= nEB + 1		# Number of nodes in width / y direction
	
	''' Beam shape '''
	shape_upchord 	= 1				# Shape of the upper chord
	shape_lowchord 	= 1				# Shape of the lower chord
									# "1" = straight line with / without slope
									# "2" = two straight lines with a kink
									# "3" = straight line followed by an circular arc segment followed by a straight line
	''' Upper chord ''' 
	#1. straight line with / without slope	
	if shape_upchord == 1:
		t_up1		= -H*0.5		# y coordinate of the upper chord at x = 0
		m_up1 		= 0.			# Slope of the upper chord in [-]
		t_up2 		= t_up1 + m_up1*L	# y coordinate of the upper chord at x = L
		upchord = []
		[upchord.append(m_up1 * intL * sizeEL + t_up1) for intL in range(nEL+1) ]
	#2. two straight lines with a kink
	elif shape_upchord == 2:
		t_up1 		= -H*0.5		# y coordinate of the upper chord at x = 0
		m_up1 		= -210./10420.	# Slope of the 1st part of the upper chord in [-]
		m_up2 		= 210./10420.	# Slope of the 2nd part of the upper chord in [-]
		nEL_kink_up = int(nEL/2)	# Number of the element after which the kink occurs
		t_up2 = t_up1 + m_up1*nEL_kink_up*sizeEL + m_up2*(nEL-nEL_kink_up)*sizeEL 	# y coordinate of the upper chord at L = 0
		upchord = []
		[upchord.append(m_up1 * intL * sizeEL + t_up1) for intL in range(int(nEL_kink_up+1)) ]
		[upchord.append(m_up1 * nEL_kink_up * sizeEL + t_up1 + (m_up2 * intL * sizeEL)) for intL in range(1, int(nEL+1-nEL_kink_up)) ]	
	#3. straight line followed by an circular arc segment followed by a straight line
	elif shape_upchord == 3:
		t_up1 		= -H*0.5		# y coordinate of the upper chord at x = 0
		t_up2 		= -H*0.5		# y coordinate of the upper chord at x = L
		m_up1 		= -1989./6936.3	# Slope of the 1st part of the upper chord in [-]
		m_up2 		= 1989./6936.3	# Slope of the 3rd part of the upper chord in [-]
		nEL_kink_up1 = 0.39*nEL		# Number of the element after which the 1st kink occurs
		nEL_kink_up2 = 0.61*nEL		# Number of the element after which the 2nd kink occurs
		radius_arc_o = 12940.3		# Radius of the arc in [mm]
		fish_up = 0					# Fish-bellied girder: 0 = no, 1 = yes
		x_MP_arc_o 	= L/2			# x coordinate of the center of the arc
		y_MP_arc_o 	= 10328.1		# y coordinate of the center of the arc
		upchord = []
		[upchord.append(m_up1 * intL * sizeEL + t_up1) for intL in range(int(nEL_kink_up1)) ]
		if fish_up == 0:
			[upchord.append(-math.sqrt(radius_arc_o**2 - (abs(nEL_kink_up1 - x_MP_arc_o/sizeEL + intL) * sizeEL)**2) + y_MP_arc_o) for intL in range(0, int(round((nEL_kink_up2+1-nEL_kink_up1),0))) ]
		else:
			[upchord.append( math.sqrt(radius_arc_o**2 - (abs(nEL_kink_up1 - x_MP_arc_o + intL) * sizeEL)**2) + y_MP_arc_o) for intL in range(0, int(round((nEL_kink_up2+1-nEL_kink_up1),0))) ]
		[upchord.append(t_up2 - (m_up2 * (nEL-nEL_kink_up2) * sizeEL) + (m_up2 * intL * sizeEL)) for intL in range(1, int(nEL+1-nEL_kink_up2)) ]
	
	''' Lower chord ''' 
	#1. straight line with / without slope
	if shape_lowchord == 1:
		t_u1 		= H*0.5			# y coordinate of the lower chord at x = 0
		m_u1 		= 0				# Slope of the lower chord in [-]
		t_u2 		= t_u1 + m_u1*L	# y coordinate of the lower chord at x = 0
		lowchord = []
		[lowchord.append(m_u1 * intL * sizeEL + t_u1) for intL in range(nEL+1) ]
	#2. two straight lines with a kink
	elif shape_lowchord == 2:
		t_u1 		= H*0.5			# y coordinate of the upper chord at x = 0
		m_u1 		= -0.1			# Slope of the 1st part of the upper chord in [-]
		m_u2 		= 0.1			# Slope of the 2nd part of the upper chord in [-]
		nEL_kink_low = int(nEL/2)	# Number of the element after which the kink occurs
		t_u2 = t_u1 + m_u1*nEL_kink_low*sizeEL + m_u2*(nEL-nEL_kink_low)*sizeEL 	# y coordinate of the upper chord at L = 0
		lowchord = []
		[lowchord.append(m_u1 * intL * sizeEL + t_u1) for intL in range(nEL_kink_low+1) ]
		[lowchord.append(m_u1 * nEL_kink_low * sizeEL + t_u1 + (m_u2 * intL * sizeEL)) for intL in range(1, int(nEL+1-nEL_kink_low)) ]
	#3. straight line followed by an circular arc segment followed by a straight line
	elif shape_lowchord == 3:
		t_u1 		= H*0.5			# y coordinate of the upper chord at x = 0
		t_u2 		= 1152			# y coordinate of the lower chord at x = 0
		m_u1 		= 828.8/11915.8	# Slope of the 1st part of the upper chord in [-]
		m_u2 		= -337.8/12089	# Slope of the 2nd part of the upper chord in [-]
		nEL_kink_low1 = nEL*0.45	# Number of the element after which the 1st kink occurs
		nEL_kink_low2 = nEL*0.54	# Number of the element after which the 2nd kink occurs
		radius_arc_low = 26829		# Radius of the arc in [mm]
		fish_low = 1				# Fish-bellied girder: 0 = no, 1 = yes
		x_MP_arc_u 	= 13724.2		# x coordinate of the center of the arc
		y_MP_arc_u 	= -25330.5		# y coordinate of the center of the arc
		lowchord = []
		[lowchord.append(m_u1 * intL * sizeEL + t_u1) for intL in range(int(nEL_kink_low1)) ]
		if fish_low == 0:
			[lowchord.append(-math.sqrt(radius_arc_low**2 - (abs(nEL_kink_low1 - x_MP_arc_u/sizeEL + intL) * sizeEL)**2) + y_MP_arc_u) for intL in range(0, int(round((nEL_kink_low2+1-nEL_kink_low1),0))) ]
		else:
			[lowchord.append( math.sqrt(radius_arc_low**2 - (abs(nEL_kink_low1 - x_MP_arc_u/sizeEL + intL) * sizeEL)**2) + y_MP_arc_u) for intL in range(0, int(round((nEL_kink_low2+1-nEL_kink_low1),0))) ]
		[lowchord.append(t_u2 - (m_u2 * (nEL-nEL_kink_low2) * sizeEL) + (m_u2 * intL * sizeEL)) for intL in range(1, int(nEL+1-nEL_kink_low2)) ]	
	
	''' Lateral bracing '''			# a maximum number of 10 supporting points is possible
	Brace_type 		= 'none'		# none, NodestoGround oder springs+beam
	if Brace_type !=	'none':
		pos_brace 	 	= [0.33,0.67]				# Position of supporting points relative to the lenght L = 1
		min_nNH_brace 	= 1							# Starting from which node from the upper chord are the springs applied, minimum 1, maximum nEH
		max_nNH_brace 	= 2							# Up to which node from the upper chord are the springs applied, minimum 2, maximum nEH+1

		# NodestoGround
		if Brace_type ==	'NodestoGround':
			K_NodestoGround = [11089.,6238.,11089.]	# List of spring stiffnesses in y direction of the supports of the horizontal bracing
		
		# Springs+beam 								# The bracing is represented by an equivalent rectangular beam with width B_beam_brace in y direction, height H_beam_brace in z direction and span = L
		elif Brace_type == 'springs+beam':			# Springs connecting the beam with the equivalent beam of the bracing
			K_beam_brace = [18300.,18300.]			# List of spring stiffnesses in y direction of the purlins, coupling of the beam with the equivalent beam representing the bracing
			K_epurlin_0	= 21200.					# Spring stiffnesses in y direction of the eaves purlin at x = 0, coupling of the beam with the equivalent beam representing the bracing
			K_epurlin_L	= 21200.					# Spring stiffnesses in y direction of the eaves purlin at x = L, coupling of the beam with the equivalent beam representing the bracing
			B_beam_brace = 3000.					# Width of the equivalent rectangular beam of the bracing in [mm]
													# Moment of inertia I_z is computed with H_beam_brace = 100 mm
			E_beam_brace = 10000.					# Modulus of elasticity of the equivalent rectangular beam of the bracing in [N/mm�]
			dist_beam_brace	= [0, 1, 2, 1, 0]		# Distances of the i-th computed beam to the bracing, 0 = Beam is part of the bracing, +1 = +1 bay -> K_beam_brace and K_epurlin are adapet automatically

			# Modify the stiffness of the purlins depending on the distance to the bracing
			for intL in range(len(K_beam_brace)):
				if dist_beam_brace[i] == 0:
					K_beam_brace[intL] = 10**9		# Beam is part of the bracing, infinite stiffness is assumed
				else:
					K_beam_brace[intL] = K_beam_brace[intL]/dist_beam_brace[i]
			if dist_beam_brace[i] == 0:
				K_epurlin_0 = 10**9					# Beam is part of the bracing, infinite stiffness is assumed
				K_epurlin_L = 10**9					# Beam is part of the bracing, infinite stiffness is assumed
			else:
				K_epurlin_0 = K_epurlin_0/dist_beam_brace[i]
				K_epurlin_L = K_epurlin_L/dist_beam_brace[i]
	
	''' Imperfections '''
	if type_var == 'manu':
		Imp_type = 'equi'				# 'equi' = sinusoidal equivalent imperfections, 'meas' = measured imperfections, 'none'
	elif type_var == 'auto': 
		if scatter_geomImp_y == 1 or scatter_geomImp_z == 1:
			Imp_type = 'equi'
		else:
			Imp_type = 'none'
			
	Imp_shape = 'bow+twist'				# Which imperfections are applied:	bow, twist, bow+twist, none
	
	# Equivalent imperfections
	if Imp_type == 'equi' and type_var == 'auto' and (scatter_geomImp_y == 1 or scatter_geomImp_z == 1):
		if scatter_geomImp_y == 1:
			ey1_equi 	= data_Imp[i][0] * L
			ey2_equi 	= 0
			ez1_equi 	= 0
			ez2_equi 	= 0
			ethetaSupp_0_equi = data_Imp[i][2]
			ethetaSupp_L_equi = data_Imp[i][3]
			ethetaMid1_equi = data_Imp[i][1] * L / H
			ethetaMid2_equi = 0
			nsin_Imp2	= 2
		else:
			ey1_equi 	= 0
			ey2_equi 	= 0
			ez1_equi 	= data_Imp[i][0] * L
			ez2_equi 	= 0
			ethetaSupp_0_equi = 0
			ethetaSupp_L_equi = 0
			ethetaMid1_equi = 0
			ethetaMid2_equi = 0
			nsin_Imp2	= 2
	elif Imp_type == 'equi':
		ey1_equi 	= L/1000			# Global bow imperfection in y direction in [mm]
		ey2_equi 	= 0					# Local bow imperfection in y direction in [mm], + -> sinusoidal half wave closest to x = 0 in +y direction
		ez1_equi 	= 0#L/1000			# Global bow imperfection in z direction in [mm]
		ez2_equi 	= 0					# Local bow imperfection in z direction in [mm], + -> sinusoidal half wave closest to x = 0 in +y direction
		ethetaSupp_0_equi = 0#1./200.	# Twist imperfection at Support 0 in [-], + -> upper edge in +y direction
		ethetaSupp_L_equi = 0#1./200.	# Twist imperfection at Support L in [-], + -> upper edge in +y direction
		ethetaMid1_equi = L/(1500*H)	# Additonal global twist imperfection at midspan in [-], + -> upper edge in +y direction
		ethetaMid2_equi = 0				# Local twist imperfection in [-], + -> upper edge of sinusoidal half wave closest to x = 0 in +y direction
		nsin_Imp2	= 2					# Number of sinusoidal half waves of local imperfections, do NOT set equal 0	
	
	if Imp_type == 'equi':
		# Imperfections of edge nodes
		ey = [0] * nNL
		ez = [0] * nNL
		etheta = [0] * nNL
		for intL in range(nNL):
			ey[intL] =  math.sin(math.pi*intL*sizeEL/L)*ey1_equi + math.sin(math.pi*intL*sizeEL/L*nsin_Imp2)*ey2_equi
			ez[intL] =  math.sin(math.pi*intL*sizeEL/L)*ez1_equi + math.sin(math.pi*intL*sizeEL/L*nsin_Imp2)*ez2_equi
			etheta[intL] =  math.sin(math.pi*intL*sizeEL/L)*ethetaMid1_equi + math.sin(math.pi*intL*sizeEL/L*nsin_Imp2)*ethetaMid2_equi + ethetaSupp_0_equi - (ethetaSupp_0_equi - ethetaSupp_L_equi)*intL/nEL
		# Imperfections of intermediate nodes
		ey_IN = [0] * nEL
		ez_IN = [0] * nEL
		etheta_IN = [0] * nEL
		for intL in range(nEL):
			ey_IN[intL]		=  math.sin(math.pi*(intL+0.5)*sizeEL/L)*ey1_equi + math.sin(math.pi*(intL+0.5)*sizeEL/L*nsin_Imp2)*ey2_equi
			ez_IN[intL]		=  math.sin(math.pi*(intL+0.5)*sizeEL/L)*ez1_equi + math.sin(math.pi*(intL+0.5)*sizeEL/L*nsin_Imp2)*ez2_equi
			etheta_IN[intL] =  math.sin(math.pi*(intL+0.5)*sizeEL/L)*ethetaMid1_equi + math.sin(math.pi*(intL+0.5)*sizeEL/L*nsin_Imp2)*ethetaMid2_equi + ethetaSupp_0_equi - (ethetaSupp_0_equi - ethetaSupp_L_equi)*(intL+0.5)/nEL		
	# Measured or arbitrary imperfections from .txt file
	elif Imp_type == 'meas':
		datei = open('data_Imp_meas_'+str(i)+'.txt', 'r')							# txt with measured/arbitrary imperfections is read in
		ey = [1000*float(column) for column in datei.readline().split(',')]		# Measured bow imperfections at edge nodes of C3D20R elements
		ey_IN = [1000*float(column) for column in datei.readline().split(',')]	# Measured bow imperfections at intermediate nodes of C3D20R elements
		ez = [0] * nNL
		ez_IN = [0] * nEL
		for intL in range(nNL):
			ez[intL] =  math.sin(math.pi*intL*sizeEL/L)*ez1_equi + math.sin(math.pi*intL*sizeEL/L*nsin_Imp2)*ez2_equi
			ez_IN[intL]		=  math.sin(math.pi*(intL+0.5)*sizeEL/L)*ez1_equi + math.sin(math.pi*(intL+0.5)*sizeEL/L*nsin_Imp2)*ez2_equi
		etheta = [float(column) for column in datei.readline().split(',')]		# Measured twist imperfections at edge nodes of C3D20R elements	
		etheta_IN = [float(column) for column in datei.readline().split(',')]	# Measured twist imperfections at intermediate nodes of C3D20R elements	
		del column
	
	########################################
	''' Boundary Conditions & Loading '''
	########################################
	''' Supports '''
	Supp_type 		= 'ideal'		# 'ideal' = ideal fork support at the cross-section centroid, 'real' = Fork support with vertical support at the lower edge of the beam with length nL_supp_L
	nEL_supp_el_0   = 2				# Length of the support zone close to x = 0 where compression plasticising and compression and tensile failure are neglected, in number of elements
	nEL_supp_el_L   = nEL_supp_el_0	# Length of the support zone close to x = L where compression plasticising and compression and tensile failure are neglected, in number of elements
	nEL_supp_shear_0= 0				# Length of the support zone close to x = 0 where shear failure is neglected, in number of elements
	nEL_supp_shear_L= nEL_supp_shear_0	# Length of the support zone close to x = L where shear failure is neglected, in number of elements
	
	# For Supp_type = 'real'
	nEL_supp_0		= 2				# Lenght of support area close to x = 0 in number of elements
	nEL_supp_L		= 2				# Lenght of support area close to x = L in number of elements
	Pos_supp_0		= 150.			# x position of the support close to x = 0 in [mm]
	Pos_supp_L		= float(L-150)	# x position of the support close to x = L in [mm]
	
	''' Loading '''					# Axial actions, Torsional actions and Actions causing My can be applied to the model in a maximum of 3 steps - can be applied in the same step as well
									# Some actions cannot be combined: Nx and My_const can only be combined if both are force controlled
	# Axial actions
	Step_n_Nx		= 0				# 0 = step is not created, 1 = 1st load step, 2 = 2nd load step , 3rd load step
	Load_type_Nx	= 'w'			# 'F' = force controled, 'Exp_LTB' = force controled only on a height of 2 elements, 'w' = displacement controlled application of axial actions
	Nx 				= 1700000.		# Normal force in [N], compression positive, applied as a surface traction over the whole end surfaces (no shear rigidity)
	wx				= 6.			# Applied displacement at x = L  in [mm], compression positive, applied as a deformation over the whole end surfaces which induces a full shear rigidity
	ey_Nx_0			= 0.			# Eccentricity to cross section centroid in y direction at x = 0 in [mm], negative values create bending analogous to positive imperfections
	ey_Nx_L			= 0.			# Eccentricity to cross section centroid in y direction at x = L in [mm], negative values create bending analogous to positive imperfections
	ez_Nx_0			= 0.			# Eccentricity to cross section centroid in z direction at x = 0 in [mm], negative values create bending analogous to positive imperfections
	ez_Nx_L			= 0.			# Eccentricity to cross section centroid in z direction at x = L in [mm], negative values create bending analogous to positive imperfections
	M_mu_supp		= 0.			# Frictional moment of the joint of the axial force at x = 0 and L, around y and z, in [Nmm]
	
	# Torsional actions
	Step_n_Mx		= 0				# 0 = step is not created, 1 = 1st load step, 2 = 2nd load step , ...
	Load_type_Mx	= 'w'			# 'F' = force controled, 'w' = displacement controlled application of torsinal actions
	Mx				= 1000000.		# Torsional moment in [Nmm], at location nEL_Mx
	rotx			= 0.05			# Applied displacement in [rad] at location nEL_Mx
	Pos_Mx			= L/2			# x position of the Torsional action
	
	# Actions causing My
	Step_n_My		= 1				# 0 = step is not created, 1 = 1st load step, 2 = 2nd load step , ...
	Load_type_My	= 'w'			# 'F' = force controled, 'w' = displacement controlled application of torsinal actions
	Moment			= 'My_const'	# Only one of the following can be applied:		Qz_3P 		= Point load at midspan - 3-point bending
									#												Qz_3Ps 		= Point load as surface traction at midspan - 3-point bending - only possible with Load_type_My = 'F'
									# 												Qz_3P_area 	= Point load at midspan distributed over an area, width = B, length = nEL_area - 3-point bending
									#												Qz_4P 		= 2 equal point loads in the third points - 4-point bending
									#												Qz_4Ps 		= 2 equal point loads as surface traction in the third points - 4-point bending	- only possible with Load_type_My = 'F'							
									#												Qz_4P_area 	= 2 equal point loads in the third points distributed over an area, width = B, length = L_area - 4-point bending
									# 												Qz_4P_area_T= 2 equal point loads in the third points distributed over an area aplied with a traverse, width = B, length = L_area - 4-point bending
									# 												My_const	= End moments causing constant My END = Endmomente
									# 												qz 			= Constant line load as surface traction at upper beam edge - only possible with Load_type_My = 'F' - refering to the actual length of the upper chord
									#												Exp_LTB		= Same as Qz_3P_area + additional horizontal force to account for vertical rotation of the loading cylinder in the experiments	
	F 				= 0.			# Force or moment in [N] for Qz_3P & Qz_4P or [Nmm] for My_const or [N/mm] for qz - for Qz_4P, each of the two loads = F, for My_const a negative value induces positive bending moment My
	w				= -0.04			# Displacement or rotation in [mm] or [rad], for My_const a negative value induces positive bending moment My, cannot be apllied for qz
	# For Qz						# Does not apply for surface tractions
									# Only the reference points are moved, for Qz_area the load application area is always at the upper chord
	nEL_area		= 2				# Lenght of load area in number of elements for Qz_3P_area
	ey_Qz			= 0.			# Eccentricity to cross section centroid in y direction in [mm], additonally to imperfections, not valid for area loading (Qz_3P_area & Qz_4P_area) when Load_type_My = 'F'
	ez_Qz			= -H/2			# Eccentricity to cross section centroid in z direction in [mm], additonally to imperfections, -H/2 = upper edge of the beam
	Pos_Qz_3P		= L/2			# x position of the action Pos_Qz_3P, in [mm]
	Pos_Qz_4P_0		= L/3			# x position of the action of Pos_Qz_4P closer to x = 0, in [mm]
	Pos_Qz_4P_L		= 2*L/3			# x position of the action of Pos_Qz_4P closer to x = L, in [mm]
	# For Qz_3P_area
	L_screws_load	= 0.			# Lenght of screws at load introduction, in this area tensile failure is neglected; vertical loads are induced over this height
	nEL_screws_load  = 2			# Lenght of loading area where tensile failure is neglected, in number of elements
	# For Exp_LTB
	ey_Exp_LTB		= 0.			# Distance of the upper hinge of the cylinder for vertical loading to its lower hinge in y direction in [mm]
	ez_Exp_LTB		= 0.			# Distance of the upper hinge of the cylinder for vertical loading to its lower hinge in z direction in [mm]
	
	########################################
	''' Material'''
	########################################
	elem_deletion	= 0				# With element deletion when tensile strengths in KaReMo are exceeded, 0 = no, 1 = yes
	
	''' Basic material '''
	Material		= 'Timber'		# Timber, Concrete, Steel - for Steel or Concrete material parameters have to be modified in line 407 and 412		
	Material_grade	= 'user'		# GL24h, GL28c/h, GL30c/h, GL32c/h, GL75, user - new timber materials can be added in the material table in line 278
	Fraktile		= '5'			# 5 = 5 % fractile and mean are implemented for timber. mean is only applied for stiffnesses E and G. Strenghts are always characteristic values, as the codes don't give any values here.
	Plasticising	= 'elpl'		# el = elastic or elpl = elastoplastic are implemented for timber. For elpl von Mises is used, which is usually alright for beams as relevant stresses and strains only occur in one material directions. 
	f_c0pl			= 0.75			# Compression sress in 1 direction when plasticising starts / proportionality limit = f_c0pl * f_c0
	e_c0pl			= 0.40			# Plastic strain when the compression strength f_c0 is reached = e_c0pl * f_c0 / E_1
	E_1pl			= 50.			# For UMAT: Modulus of elasticity for softening: E_1ges = (1/E_1+1/E_1pl)^-1
	f_v12pl1		= 0.30			# For UMAT: Shear stress in 12 direction when plasticising starts / proportionality limit = f_vpl1 * f_v
	f_v12pl2		= 0.30			# For UMAT: Shear stress in 12 direction when softening starts / proportionality limit = f_vpl2 * f_v
	e_v21pl			= 0.00			# For UMAT: Plastic strain when the shear strength f_v is reached = e_vpl * f_v / G_12	
	G_12pl			= 975.			# For UMAT: Shear modulus for softening: G_12ges = (1/G_12+1/G_12pl)^-1		
	f_v13pl1		= 0.30			# For UMAT: Shear stress in 12 direction when plasticising starts / proportionality limit = f_vpl1 * f_v
	f_v13pl2		= 0.30			# For UMAT: Shear stress in 12 direction when softening starts / proportionality limit = f_vpl2 * f_v
	e_v31pl			= 0.00			# For UMAT: Plastic strain when the shear strength f_v is reached = e_vpl * f_v / G_12	
	G_13pl			= 975.			# For UMAT: Shear modulus for softening: G_12ges = (1/G_12+1/G_12pl)^-1		
	f_v23pl1		= 0.55			# For UMAT: Shear stress in 12 direction when plasticising starts / proportionality limit = f_vpl1 * f_v
	f_v23pl2		= 0.55			# For UMAT: Shear stress in 12 direction when softening starts / proportionality limit = f_vpl2 * f_v
	e_v23pl			= 0.00			# For UMAT: Plastic strain when the shear strength f_v is reached = e_vpl * f_v / G_12	
	G_23pl			= 50.			# For UMAT: Shear modulus for softening: G_12ges = (1/G_12+1/G_12pl)^-1
	Material_law	= 'UMAT'		# Plasticising with 'UMAT' or 'VM' = von Mises. If UMAT is used for Edge lamellae UMAT_tedge.for and UMAT_cedge.for is needed
	UMAT_format 	= 'dll'			# File format of UMAT. Either 'for', 'f' or 'dll'. dll can only be used if no edge lamellaes are defined, as no name variation is possible to account for different UMATs
									# The UMAT file must be named either: 'UMAT.for', 'UMAT.f' or 'standardU.dll'
									# A 'dll' can be run without a Fortran compiler	
									
	''' Optional - Edge lamellae with different material '''	# Cannot be combined with steel and concrete
	nEdgelam_Tens	= 0				# Number of elements of the tension edge lamellae in height
	MatEdgeLam_Tens	= 'GL75'		# Material grade of the tension edge lamellaes, fractile at plasticising are chosen in the same way as for the basic material
	nEdgelam_Comp	= 0				# Number of elements of the compression edge lamellae in height
	MatEdgeLam_Comp	= 'GL75'		# Material grade of the compression edge lamellaes, fractile at plasticising are chosen in the same way as for the basic material
	f_c0pl_edge		= 0.75			# Compression sress in 1 direction when plasticising starts / proportionality limit = f_c0pl * f_c0
	e_c0pl_edge		= 0.40			# Plastic strain when the compression strength f_c0 is reached = e_c0pl * f_c0 / E_1
	E_1pl_edge		= 50.			# For UMAT: Modulus of elasticity for softening: E_1ges = (1/E_1+1/E_1pl)^-1
	f_v12pl1_edge	= 0.60			# For UMAT: Shear stress in 12 direction when plasticising starts / proportionality limit = f_vpl1 * f_v
	f_v12pl2_edge	= 0.70			# For UMAT: Shear stress in 12 direction when softening starts / proportionality limit = f_vpl2 * f_v
	e_v21pl_edge	= 0.30			# For UMAT: Plastic strain when the shear strength f_v is reached = e_vpl * f_v / G_12	
	G_12pl_edge		= 1.			# For UMAT: Shear modulus for softening: G_12ges = (1/G_12+1/G_12pl)^-1
	f_v13pl1_edge	= 0.60			# For UMAT: Shear stress in 12 direction when plasticising starts / proportionality limit = f_vpl1 * f_v
	f_v13pl2_edge	= 0.70			# For UMAT: Shear stress in 12 direction when softening starts / proportionality limit = f_vpl2 * f_v
	e_v31pl_edge	= 0.30			# For UMAT: Plastic strain when the shear strength f_v is reached = e_vpl * f_v / G_12	
	G_13pl_edge		= 1.			# For UMAT: Shear modulus for softening: G_12ges = (1/G_12+1/G_12pl)^-1
	f_v23pl1_edge	= 0.60			# For UMAT: Shear stress in 12 direction when plasticising starts / proportionality limit = f_vpl1 * f_v
	f_v23pl2_edge	= 0.70			# For UMAT: Shear stress in 12 direction when softening starts / proportionality limit = f_vpl2 * f_v
	e_v23pl_edge	= 0.30			# For UMAT: Plastic strain when the shear strength f_v is reached = e_vpl * f_v / G_12	
	G_23pl_edge		= 1.			# For UMAT: Shear modulus for softening: G_12ges = (1/G_12+1/G_12pl)^-1
	Material_law_edge= 'VM'			# Plasticising with 'UMAT' or 'VM' = von Mises. If UMAT is used for Edge lamellae UMAT_tedge.for and UMAT_cedge.for is needed
	UMAT_format_edge= 'dll'			# File format of UMAT. Either 'for' or 'f'
									# The UMAT file must be named either: 'UMAT.for' or 'UMAT.f'	
									
	''' Material orientation '''
	Mat_axis_vert	= 'R'			# Material direction perpendicular to the lower edge of the beam, 'R' or 'T'
	
	########################################
	''' Analysis'''
	########################################
	''' Analysis type '''
	Cal_type		= 'GMNIA'		# GMNIA = Geometrically and materially nonlinear abalysis with imperfections, LBA = linear buckling analysis, the reference total load is 1 N or Nmm
									# For LBA all loads have to be applied in step one with their total value ranging from 0 to 1
	''' Solver settings '''			# Only the steps defined in are created
	startInc_One 	= 0.05			# Increment size at the beginning of step one
	minInc_One		= 0.001			# Minimum increment size in step one
	maxInc_One		= 0.05			# Maximum increment size in step one
	maxnInc_One		= 1000			# Maxmimum number of increments in step one
	startInc_Two 	= 0.05			# Increment Size at the beginning of Step two
	minInc_Two		= 0.001			# Minimum increment size in step two
	maxInc_Two		= 0.05			# Maximum increment size in step two
	maxnInc_Two 	= 1000			# Maxmimum number of increments in step two
	startInc_Three 	= 0.05			# Increment Size at the beginning of Step three
	minInc_Three	= 0.001			# Minimum increment size in step three
	maxInc_Three	= 0.05			# Maximum increment size in step three
	maxnInc_Three 	= 1000			# Maxmimum number of increments in step three
	Sol_tech 		= QUASI_NEWTON	# FULL_NEWTON (default), QUASI_NEWTON (for convergence problems with nonlinear materials)
	
	########################################
	''' Saving the model '''
	########################################
	strPath = "G:\Abaqus\work_directory"	# Set storage directory
	modelName = str(int(round(L,0)))+'x'+str(int(round(H,0)))+'x'+str(int(round(B,0)))+'_'+Material_grade+'-'+Fraktile+'-'+Plasticising+'_'+str(Cal_type)		# Model name
	modelName_save = str(int(round(L,0)))+'x'+str(int(round(H,0)))+'x'+str(int(round(B,0)))+'_'+Material_grade+'-'+Fraktile+'-'+Plasticising+'_'+str(Cal_type)+'_'+str(i)		# Model name
	myModel = mdb.Model(modelType=STANDARD_EXPLICIT, name=modelName)	# Model type
	
	########################################
	''' Material data '''
	########################################
	''' Timber '''
	if Material == 'Timber':
		# Orthotropic material with the three directions:	1 = Longitudinal	in x direction
		#													2 = Radial			in y direction
		#													3 = Tangential		in z direction
		
		# Material database for timber 
		# 						'Name'			E_1		E_2		E_3		v_12	v_13	v_23	G_23	G_12	G_13	f_c0	f_t0	f_v12	f_v13	f_v23	
		Materialdb = []
		[Materialdb.append((	'user',			11500,	300,	300,	0.50,	0.50,	0.50,	100,	650,	650,	24.0,	24.0,	3.5,	3.5,	1.2))]		# can be modified 
		[Materialdb.append((	'GL24h_mean',	11500,	300,	300,	0.50,	0.50,	0.50,	100,	650,	650,	40.0,	33.0,	5.3,	5.3,	1.6))]		# according to Diss_Töpler
		[Materialdb.append((	'GL24h_5',		9600,	250,	250,	0.50,	0.50,	0.50,	54,		540,	540,	24.0,	24.0,	3.5,	3.5,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL24h_95',		13900,	360,	360,	0.50,	0.50,	0.50,	79,		790,	790,	43.0,	40.0,	9.0,	9.0,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL24c_mean',	11000,	300,	300,	0.50,	0.50,	0.50,	65,		650,	650,	21.5,	24.0,	3.5,	3.5,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL24c_5',		9100,	250,	250,	0.50,	0.50,	0.50,	54,		540,	540,	21.5,	24.0,	3.5,	3.5,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL28h_mean',	12600,	300,	300,	0.50,	0.50,	0.50,	65,		650,	650,	28.0,	28.0,	3.5,	3.5,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL28h_5',		10500,	250,	250,	0.50,	0.50,	0.50,	54,		540,	540,	28.0,	28.0,	3.5,	3.5,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL28c_mean',	12500,	300,	300,	0.50,	0.50,	0.50,	65,		650,	650,	24.0,	28.0,	3.5,	3.5,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL28c_5',		10400,	250,	250,	0.50,	0.50,	0.50,	54,		540,	540,	24.0,	28.0,	3.5,	3.5,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL30c_mean',	13000,	300,	300,	0.50,	0.50,	0.50,	100,	650,	650,	45.0,	41.0,	5.3,	5.3,	1.6))]		# according to Diss_Töpler
		[Materialdb.append((	'GL30c_5',		10800,	250,	250,	0.50,	0.50,	0.50,	54,		540,	540,	24.5,	30.0,	3.5,	3.5,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL30h_mean',	13600,	300,	300,	0.50,	0.50,	0.50,	65,		650,	650,	30.0,	30.0,	3.5,	3.5,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL30h_5',		11300,	250,	250,	0.50,	0.50,	0.50,	54,		540,	540,	30.0,	30.0,	3.5,	3.5,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL32c_mean',	13500,	300,	300,	0.50,	0.50,	0.50,	65,		650,	650,	24.5,	30.0,	3.5,	3.5,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL32c_5',		11200,	250,	250,	0.50,	0.50,	0.50,	54,		540,	540,	24.5,	30.0,	3.5,	3.5,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL32h_mean',	14200,	300,	300,	0.50,	0.50,	0.50,	65,		650,	650,	30.0,	32.0,	3.5,	3.5,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL32h_5',		11800,	250,	250,	0.50,	0.50,	0.50,	54,		540,	540,	30.0,	32.0,	3.5,	3.5,	1.2))]		# E, G and f according to DIN EN 14080:2013, v means calculated with values of Neuhaus (2017) for u = 12 %
		[Materialdb.append((	'GL75_mean',	16800,	900,	900,	0.30,	0.50,	0.20,	150,	900,	1000,	77.0,	90.0,	100,	100,	100))]		# according to Diss_Töpler
		[Materialdb.append((	'GL75_5',		15300,	400,	400,	0.30,	0.50,	0.20,	76,		760,	760,	59.4,	75.0,	4.5,	4.5,	3.8))]		# E, G and f according to ETA-14/0354 (2018), v according to Töpler & Kuhlmann (2022) for u = 5.5 %
		
		# No scattering material values
		if type_var == 'manu' or (KaReMo == 0 and KaReMo_v == 0 and KaReMo_y == 0):
		
			# Read Material from database
			for intZ in range(1,4):
				if intZ == 1:							# Basic material
					if Material_grade == 'user':
						myMaterial = Material_grade
					else:
						myMaterial = Material_grade+'_'+Fraktile
				elif intZ == 2 and nEdgelam_Tens != 0:	# Tension edge lamellae
					if MatEdgeLam_Tens == 'user':
						myMaterial = MatEdgeLam_Tens
					else:
						myMaterial = MatEdgeLam_Tens+'_'+Fraktile
				elif intZ == 3 and nEdgelam_Comp != 0:	# Compression edge lamellae
					if MatEdgeLam_Tens == 'user':
						myMaterial = MatEdgeLam_Comp
					else:
						myMaterial = MatEdgeLam_Comp+'_'+Fraktile		
					
				for intM in range(len(Materialdb)):
					row = Materialdb[intM]
					if row[0] == myMaterial:
						break
				
				# Assign material values
				E_1		= row[1]			#[N/mm^2]
				E_2		= row[2]			#[N/mm^2]
				E_3		= row[3]			#[N/mm^2]
				v_12	= row[4]			#[-]
				v_13	= row[5]			#[-]
				v_23	= row[6]			#[-]			
				G_23	= row[7]			#[N/mm^2]
				G_12	= row[8]			#[N/mm^2]
				G_13	= row[9]			#[N/mm^2]
				f_c0	= row[10]			#[N/mm^2]
				f_t0	= row[11]			#[N/mm^2]
				f_v12	= row[12]			#[N/mm^2
				f_v13	= row[13]			#[N/mm^2]
				f_v23	= row[14]			#[N/mm^2]
				v_21	= (v_12*E_2)/E_1
				v_31	= (v_13*E_3)/E_1
				v_32	= (v_23*E_3)/E_2
				Delta	= 1/(1-v_12*v_21-v_13*v_31-v_23*v_32-2*v_13*v_21*v_32)

				# Calculate stiffness matrix for direct input
				D1111 = round((1-v_32*v_23)*Delta*E_1)
				D1122 = round((v_13*v_32+v_12)*Delta*E_2)
				D2222 = round((1-v_13*v_31)*Delta*E_2)
				D1133 = round((v_12*v_23+v_13)*Delta*E_3)
				D2233 = round((v_21*v_13+v_23)*Delta*E_3)
				D3333 = round((1-v_12*v_21)*Delta*E_3)
				D1112 = 0.0
				D2212 = 0.0
				D3312 = 0.0
				D1212 = round(G_12)
				D1113 = 0.0
				D2213 = 0.0
				D3313 = 0.0
				D1213 = 0.0
				D1313 = round(G_13)
				D1123 = 0.0
				D2223 = 0.0
				D3323 = 0.0
				D1223 = 0.0
				D1323 = 0.0
				D2323 = round(G_23)
				
				# Define material models in Abaqus	
				# Basic material
				if intZ == 1:
					# Timber in Tension
					if Plasticising == 'el' or (Plasticising == 'elpl' and Material_law == 'VM'):
						myMaterialGL = myModel.Material(name='Timber_tens')
						myMaterialGL.Elastic(type=ANISOTROPIC,table=((D1111, D1122, D2222, D1133, D2233, D3333, D1112, D2212, D3312, D1212, D1113, D2213, D3313, D1213, D1313, D1123, D2223, D3323, D1223, D1323, D2323), ))
						#if Plasticising == 'elpl' and Material_law == 'VM':	# If edge lamellaes with different material are used, this is applied only to the edge lamellaes
						#	myMaterialGL.Plastic(table=((f_t0, 0.0), (0.01, 5*float(f_t0)/E_1)))	# This is a very rough softening criterion which sometimes allows for a bit of post failure modelling.
																									# Attention: Stresses are computed and evaluated at the integration point (not equal to the element nodes).
																									#			 Stresses at element nodes are extrapolated and might therefore violate this condition (and exceed f_t0)
																									# If this is really relevant, it is strongly recommended to implement something more accurate.
																									# The author often only applies failure criteria in post processing (see Output) because there, 
																									# values at element nodes can be analysed which usually govern in design verification.
					elif Plasticising == 'elpl' and Material_law == 'UMAT':
						myMaterialGL = myModel.Material(name='Timber_tens')
						myMaterialGL.UserMaterial(type=MECHANICAL, mechanicalConstants=(E_1, E_2, E_3, v_12, v_13, v_23, G_12, G_13, G_23, f_c0pl, -f_c0, e_c0pl, E_1pl, f_v12pl1, f_v12pl2*f_v12, e_v21pl, G_12pl, f_v13pl1, f_v13pl2*f_v13, e_v31pl, G_13pl, f_v23pl1, f_v23pl2*f_v23, e_v23pl, G_23pl, f_v12, f_t0))
						myMaterialGL.Depvar(n=62)													# Number of state variables (STATEV) defined in the UMAT
						
					# Timber in Compression
					if Plasticising == 'el' or (Plasticising == 'elpl' and Material_law == 'VM'):
						myMaterialGL = myModel.Material(name='Timber_comp')
						myMaterialGL.Elastic(type=ANISOTROPIC,table=((D1111, D1122, D2222, D1133, D2233, D3333, D1112, D2212, D3312, D1212, D1113, D2213, D3313, D1213, D1313, D1123, D2223, D3323, D1223, D1323, D2323), ))
						if Plasticising == 'elpl' and Material_law == 'VM':
							myMaterialGL.Plastic(table=((((f_c0pl+(1-f_c0pl)*sin(0.0*pi/2)))*f_c0, e_c0pl*(1-cos(0.0*pi/2))*f_c0/E_1), 
														(((f_c0pl+(1-f_c0pl)*sin(0.1*pi/2)))*f_c0, e_c0pl*(1-cos(0.1*pi/2))*f_c0/E_1), 
														(((f_c0pl+(1-f_c0pl)*sin(0.2*pi/2)))*f_c0, e_c0pl*(1-cos(0.2*pi/2))*f_c0/E_1), 
														(((f_c0pl+(1-f_c0pl)*sin(0.3*pi/2)))*f_c0, e_c0pl*(1-cos(0.3*pi/2))*f_c0/E_1), 
														(((f_c0pl+(1-f_c0pl)*sin(0.4*pi/2)))*f_c0, e_c0pl*(1-cos(0.4*pi/2))*f_c0/E_1), 
														(((f_c0pl+(1-f_c0pl)*sin(0.5*pi/2)))*f_c0, e_c0pl*(1-cos(0.5*pi/2))*f_c0/E_1), 
														(((f_c0pl+(1-f_c0pl)*sin(0.6*pi/2)))*f_c0, e_c0pl*(1-cos(0.6*pi/2))*f_c0/E_1), 
														(((f_c0pl+(1-f_c0pl)*sin(0.7*pi/2)))*f_c0, e_c0pl*(1-cos(0.7*pi/2))*f_c0/E_1), 
														(((f_c0pl+(1-f_c0pl)*sin(0.8*pi/2)))*f_c0, e_c0pl*(1-cos(0.8*pi/2))*f_c0/E_1), 
														(((f_c0pl+(1-f_c0pl)*sin(0.9*pi/2)))*f_c0, e_c0pl*(1-cos(0.9*pi/2))*f_c0/E_1), 
														(((f_c0pl+(1-f_c0pl)*sin(1.0*pi/2)))*f_c0, e_c0pl*(1-cos(1.0*pi/2))*f_c0/E_1), 
														(1.001*f_c0, 0.25)))
					elif Plasticising == 'elpl' and Material_law == 'UMAT':
						myMaterialGL = myModel.Material(name='Timber_comp')
						myMaterialGL.UserMaterial(type=MECHANICAL, mechanicalConstants=(E_1, E_2, E_3, v_12, v_13, v_23, G_12, G_13, G_23, f_c0pl, -f_c0, e_c0pl, E_1pl, f_v12pl1, f_v12pl2*f_v12, e_v21pl, G_12pl, f_v13pl1, f_v13pl2*f_v13, e_v31pl, G_13pl, f_v23pl1, f_v23pl2*f_v23, e_v23pl, G_23pl, f_v12, f_t0))
						myMaterialGL.Depvar(n=62)													# Number of state variables (STATEV) defined in the UMAT

					# Timber at load introduction - tensile failure is neglected
					if nEL_screws_load != 0 and L_screws_load != 0:
						if Plasticising == 'el' or (Plasticising == 'elpl' and Material_law == 'VM'):
							myMaterialGL = myModel.Material(name='Timber_load')
							myMaterialGL.Elastic(type=ANISOTROPIC,table=((D1111, D1122, D2222, D1133, D2233, D3333, D1112, D2212, D3312, D1212, D1113, D2213, D3313, D1213, D1313, D1123, D2223, D3323, D1223, D1323, D2323), ))
						elif Plasticising == 'elpl' and Material_law == 'UMAT':
							myMaterialGL = myModel.Material(name='Timber_load')
							myMaterialGL.UserMaterial(type=MECHANICAL, mechanicalConstants=(E_1, E_2, E_3, v_12, v_13, v_23, G_12, G_13, G_23, f_c0pl, -f_c0, e_c0pl, E_1pl, f_v12pl1, f_v12pl2*f_v12, e_v21pl, G_12pl, f_v13pl1, f_v13pl2*f_v13, e_v31pl, G_13pl, f_v23pl1, f_v23pl2*f_v23, e_v23pl, G_23pl, f_v12, 1000))
							myMaterialGL.Depvar(n=62)													# Number of state variables (STATEV) defined in the UMAT
					
					# Timber at supports - compression plasticising, compression and tensile failure are neglected
					if nEL_supp_el_0 != 0 or nEL_supp_el_L != 0:
						if Plasticising == 'el' or (Plasticising == 'elpl' and Material_law == 'VM'):
							myMaterialGL = myModel.Material(name='Timber_supp')
							myMaterialGL.Elastic(type=ANISOTROPIC,table=((D1111, D1122, D2222, D1133, D2233, D3333, D1112, D2212, D3312, D1212, D1113, D2213, D3313, D1213, D1313, D1123, D2223, D3323, D1223, D1323, D2323), ))
						elif Plasticising == 'elpl' and Material_law == 'UMAT':
							myMaterialGL = myModel.Material(name='Timber_supp')
							myMaterialGL.UserMaterial(type=MECHANICAL, mechanicalConstants=(E_1, E_2, E_3, v_12, v_13, v_23, G_12, G_13, G_23, f_c0pl, -1000, e_c0pl, E_1pl, f_v12pl1, f_v12pl2*f_v12, e_v21pl, G_12pl, f_v13pl1, f_v13pl2*f_v13, e_v31pl, G_13pl, f_v23pl1, f_v23pl2*f_v23, e_v23pl, G_23pl, f_v12, 1000))
							myMaterialGL.Depvar(n=62)													# Number of state variables (STATEV) defined in the UMAT
						
				# Tension edge lamellae
				elif intZ == 2 and nEdgelam_Tens != 0:  
					if Plasticising == 'el' or Material_law_edge == 'VM':
						myMaterialGL = myModel.Material(name='Timber_Edgelam_Tens')
						myMaterialGL.Elastic(type=ANISOTROPIC,table=((D1111, D1122, D2222, D1133, D2233, D3333, D1112, D2212, D3312, D1212, D1113, D2213, D3313, D1213, D1313, D1123, D2223, D3323, D1223, D1323, D2323), ))
						#if Plasticising == 'elpl':
						#	myMaterialGL.Plastic(table=((f_t0, 0.0), (0.01, 5*f_t0/E_1)))			# See comments line 446
					elif Plasticising == 'elpl' and Material_law_edge == 'UMAT':
						myMaterialGL = myModel.Material(name='Timber_Edgelam_Tens')
						myMaterialGL.UserMaterial(type=MECHANICAL, mechanicalConstants=(E_1, E_2, E_3, v_12, v_13, v_23, G_12, G_13, G_23, f_c0pl, -f_c0, e_c0pl_edge, E_1pl_edge, f_v12pl1_edge, f_v12pl2_edge*f_v12, e_v21pl_edge, G_12pl_edge, f_v13pl1_edge, f_v13pl2_edge*f_v13, e_v31pl_edge, G_13pl_edge, f_v23pl1_edge, f_v23pl2_edge*f_v23, e_v23pl_edge, G_23pl_edge, f_v12, f_t0))
						myMaterialGL.Depvar(n=62)													# Number of state variables (STATEV) defined in the UMAT
						
				# Compression edge lamellae
				elif intZ == 3 and nEdgelam_Comp != 0:
					if Plasticising == 'el' or Material_law_edge == 'VM':
						myMaterialGL = myModel.Material(name='Timber_Edgelam_Comp')
						myMaterialGL.Elastic(type=ANISOTROPIC,table=((D1111, D1122, D2222, D1133, D2233, D3333, D1112, D2212, D3312, D1212, D1113, D2213, D3313, D1213, D1313, D1123, D2223, D3323, D1223, D1323, D2323), ))
						if Plasticising == 'elpl':
							myMaterialGL.Plastic(table=((((f_c0pl_edge+(1-f_c0pl_edge)*sin(0.0*pi/2)))*f_c0, e_c0pl_edge*(1-cos(0.0*pi/2))*f_c0/E_1), 
														(((f_c0pl_edge+(1-f_c0pl_edge)*sin(0.1*pi/2)))*f_c0, e_c0pl_edge*(1-cos(0.1*pi/2))*f_c0/E_1), 
														(((f_c0pl_edge+(1-f_c0pl_edge)*sin(0.2*pi/2)))*f_c0, e_c0pl_edge*(1-cos(0.2*pi/2))*f_c0/E_1), 
														(((f_c0pl_edge+(1-f_c0pl_edge)*sin(0.3*pi/2)))*f_c0, e_c0pl_edge*(1-cos(0.3*pi/2))*f_c0/E_1), 
														(((f_c0pl_edge+(1-f_c0pl_edge)*sin(0.4*pi/2)))*f_c0, e_c0pl_edge*(1-cos(0.4*pi/2))*f_c0/E_1), 
														(((f_c0pl_edge+(1-f_c0pl_edge)*sin(0.5*pi/2)))*f_c0, e_c0pl_edge*(1-cos(0.5*pi/2))*f_c0/E_1), 
														(((f_c0pl_edge+(1-f_c0pl_edge)*sin(0.6*pi/2)))*f_c0, e_c0pl_edge*(1-cos(0.6*pi/2))*f_c0/E_1), 
														(((f_c0pl_edge+(1-f_c0pl_edge)*sin(0.7*pi/2)))*f_c0, e_c0pl_edge*(1-cos(0.7*pi/2))*f_c0/E_1), 
														(((f_c0pl_edge+(1-f_c0pl_edge)*sin(0.8*pi/2)))*f_c0, e_c0pl_edge*(1-cos(0.8*pi/2))*f_c0/E_1), 
														(((f_c0pl_edge+(1-f_c0pl_edge)*sin(0.9*pi/2)))*f_c0, e_c0pl_edge*(1-cos(0.9*pi/2))*f_c0/E_1), 
														(((f_c0pl_edge+(1-f_c0pl_edge)*sin(1.0*pi/2)))*f_c0, e_c0pl_edge*(1-cos(1.0*pi/2))*f_c0/E_1), 
														(1.001*f_c0, 0.25)))
					elif Plasticising == 'elpl' and Material_law_edge == 'UMAT':
						myMaterialGL = myModel.Material(name='Timber_Edgelam_Comp')
						myMaterialGL.UserMaterial(type=MECHANICAL, mechanicalConstants=(E_1, E_2, E_3, v_12, v_13, v_23, G_12, G_13, G_23, f_c0pl, -f_c0, e_c0pl_edge, E_1pl_edge, f_v12pl1_edge, f_v12pl2_edge*f_v12, e_v21pl_edge, G_12pl_edge, f_v13pl1_edge, f_v13pl2_edge*f_v13, e_v31pl_edge, G_13pl_edge, f_v23pl1_edge, f_v23pl2_edge*f_v23, e_v23pl_edge, G_23pl_edge, f_v12, f_t0))
						myMaterialGL.Depvar(n=62)													# Number of state variables (STATEV) defined in the UMAT
		
		# Scattering material values
		elif KaReMo == 1:
			# Read scattering material data
			if newLam == 1:
				ndata_Mat = int(H/Hlam*(L/150+54))
			else:
				ndata_Mat = int(H/Hlam*L/150)
			data_Mat = np.zeros((ndata_Mat,6))
			c = 1
			j = 0
			with open('data_beam_'+str(i)+'.txt') as f:
				for x in f:
					if c >= 2 and c < 2+ndata_Mat:
						x = x.split(',')
						if KaReMo_v == 1 and KaReMo_y == 1:
							if float(x[2]) == 0:		# No FJ - fc, ft, fv, Et_le, Et_ri, G
								data_Mat[j] = ([float(x[6]),float(x[8]),float(x[21]),float(x[26]),float(x[27]),float(x[22])])
							else:
								data_Mat[j] = ([float(x[10]),float(x[12]),float(x[21]),float(x[26]),float(x[27]),float(x[22])])
						elif KaReMo_v == 1: 
							if float(x[2]) == 0:		# No FJ - fc, ft, fv, Et_le, Et_ri, G
								data_Mat[j] = ([float(x[6]),float(x[8]),float(x[21]),float(x[7]),float(x[7]),float(x[22])])
							else:
								data_Mat[j] = ([float(x[10]),float(x[12]),float(x[21]),float(x[11]),float(x[11]),float(x[22])])
						elif KaReMo_y == 1:
							if float(x[2]) == 0:		# No FJ - fc, ft, fv, Et_le, Et_ri, G
								data_Mat[j] = ([float(x[6]),float(x[8]),float(x[21]),float(x[26]),float(x[27]),0])
							else:
								data_Mat[j] = ([float(x[10]),float(x[12]),float(x[21]),float(x[26]),float(x[27]),0])
						else:
							if float(x[2]) == 0:
								data_Mat[j] = ([float(x[6]),float(x[8]),float(x[21]),float(x[7]),float(x[7]),0])
							else:
								data_Mat[j] = ([float(x[10]),float(x[12]),float(x[21]),float(x[11]),float(x[11]),0])
						j = j + 1
					c = c + 1
			
			# Number of elements of length = 1 segment
			nELpSeg = nEL/(L/150.)				# Number of elements per segment
			if nELpSeg != round(nELpSeg,0):		# Exit if elementsize is not consistent with segment length
				error_exit = 1
				if L/150. != round(L/150.,0):
					print('Error: L is not a multiple of 150 mm')
				else:
					print('Error: nEL is not suitable')
				break
			else:
				nELpSeg = int(nELpSeg)
			
			# Pick other material values from material table
			if Material_grade == 'user':
				myMaterial = Material_grade
			else:
				myMaterial = Material_grade+'_'+Fraktile
			for intM in range(len(Materialdb)):
				row = Materialdb[intM]
				if row[0] == myMaterial:
					break
				
			# Loop over beam height
			j = 0
			k = 1
			for intH in range(int(H/Hlam)):
				# Loop over beam length
				for intL in range(int(L/150*nELpSeg)):
				# Loop over beam width
					for intW in range(2):
						# Assign material values
						if intW == 0:
							E_1	= data_Mat[j,3]		#[N/mm^2]
						else:
							E_1	= data_Mat[j,4]		#[N/mm^2]
						E_2		= row[2]			#[N/mm^2]
						E_3		= row[3]			#[N/mm^2]
						v_12	= row[4]			#[-]
						v_13	= row[5]			#[-]
						v_23	= row[6]			#[-]
						if KaReMo_v == 1:
							G_12= data_Mat[j,5]		#[N/mm^2]
							G_13= data_Mat[j,5]		#[N/mm^2]
							f_v12=data_Mat[j,2]		#[N/mm^2
							f_v13=data_Mat[j,2]		#[N/mm^2]
						else:
							G_12= row[8]			#[N/mm^2]
							G_13= row[9]			#[N/mm^2]
							f_v12=row[12]			#[N/mm^2]
							f_v13=row[13]			#[N/mm^2]						
						G_23	= row[7]			#[N/mm^2]
						f_v23	= row[14]			#[N/mm^2]
						# Timber at supports - compression plasticising, compression and tensile failure are neglected
						if intL < int(nEL_supp_el_0) or intL >= nEL-int(nEL_supp_el_L):
							f_c0= 1000				#[N/mm^2]
						else:
							f_c0= data_Mat[j,0]		#[N/mm^2]	
						# Timber at supports or load introduction - tensile failure is neglected
						if intL < int(nEL_supp_el_0) or intL >= nEL-int(nEL_supp_el_L) or (L_screws_load > 0.0 and Step_n_My != 0 and Moment == 'Qz_3P_area' and intL >= int(nEL/2)-int(nEL_screws_load/2) and intL < int(nEL/2)+int(nEL_screws_load/2) and intH < max(1,int(L_screws_load*nEH/H+1))):
							f_t0= 1000				#[N/mm^2]
						else:
							f_t0= data_Mat[j,1]		#[N/mm^2]
						v_21	= (v_12*E_2)/E_1
						v_31	= (v_13*E_3)/E_1
						v_32	= (v_23*E_3)/E_2
						Delta	= 1/(1-v_12*v_21-v_13*v_31-v_23*v_32-2*v_13*v_21*v_32)
						if intW == 1:
							if k == nELpSeg:
								j = j + 1
								k = 1
								if intL == (int(L/150*nELpSeg)-1) and newLam == 1:
									j = j + 54
							else:
								k = k + 1

						# Calculate stiffness matrix for direct input
						D1111 = round((1-v_32*v_23)*Delta*E_1)
						D1122 = round((v_13*v_32+v_12)*Delta*E_2)
						D2222 = round((1-v_13*v_31)*Delta*E_2)
						D1133 = round((v_12*v_23+v_13)*Delta*E_3)
						D2233 = round((v_21*v_13+v_23)*Delta*E_3)
						D3333 = round((1-v_12*v_21)*Delta*E_3)
						D1112 = 0.0
						D2212 = 0.0
						D3312 = 0.0
						D1212 = round(G_12)
						D1113 = 0.0
						D2213 = 0.0
						D3313 = 0.0
						D1213 = 0.0
						D1313 = round(G_13)
						D1123 = 0.0
						D2223 = 0.0
						D3323 = 0.0
						D1223 = 0.0
						D1323 = 0.0
						D2323 = round(G_23)
						
						# Define material models in Abaqus
						myMaterialGL = myModel.Material(name='Timber_H'+str(intH+1)+'_L'+str(intL+1)+'_W'+str(intW+1))
						if Material_law == 'UMAT':
							myMaterialGL.UserMaterial(type=MECHANICAL, mechanicalConstants=(E_1, E_2, E_3, v_12, v_13, v_23, G_12, G_13, G_23, f_c0pl, -f_c0, e_c0pl, E_1pl, f_v12pl1, f_v12pl2*f_v12, e_v21pl, G_12pl, f_v13pl1, f_v13pl2*f_v13, e_v31pl, G_13pl, f_v23pl1, f_v23pl2*f_v23, e_v23pl, G_23pl, f_v12, f_t0))
							if elem_deletion == 1:
								myMaterialGL.Depvar(n=62, deleteVar=9)
							else:
								myMaterialGL.Depvar(n=62)
						else:	
							myMaterialGL.Elastic(type=ANISOTROPIC,table=((D1111, D1122, D2222, D1133, D2233, D3333, D1112, D2212, D3312, D1212, D1113, D2213, D3313, D1213, D1313, D1123, D2223, D3323, D1223, D1323, D2323), ))
							if Plasticising == 'elpl' and Material_law == 'VM':
								myMaterialGL.Plastic(table=((((f_c0pl+(1-f_c0pl)*sin(0.0*pi/2)))*f_c0, e_c0pl*(1-cos(0.0*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.1*pi/2)))*f_c0, e_c0pl*(1-cos(0.1*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.2*pi/2)))*f_c0, e_c0pl*(1-cos(0.2*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.3*pi/2)))*f_c0, e_c0pl*(1-cos(0.3*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.4*pi/2)))*f_c0, e_c0pl*(1-cos(0.4*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.5*pi/2)))*f_c0, e_c0pl*(1-cos(0.5*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.6*pi/2)))*f_c0, e_c0pl*(1-cos(0.6*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.7*pi/2)))*f_c0, e_c0pl*(1-cos(0.7*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.8*pi/2)))*f_c0, e_c0pl*(1-cos(0.8*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.9*pi/2)))*f_c0, e_c0pl*(1-cos(0.9*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(1.0*pi/2)))*f_c0, e_c0pl*(1-cos(1.0*pi/2))*f_c0/E_1), 
															(f_c0+0.25*E_1pl, (f_c0/E_1+0.25))))
															
		# Scattering material values
		elif KaReMo == 0 and KaReMo_y == 1:
			# Read scattering material data
			if newLam == 1:
				ndata_Mat = int(H/Hlam*(L/150+54))
			else:
				ndata_Mat = int(H/Hlam*L/150)
			data_Mat = np.zeros((ndata_Mat,3))
			c = 1
			j = 0
			with open('data_beam_'+str(i)+'.txt') as f:
				for x in f:
					if c >= 2 and c < 2+ndata_Mat:
						x = x.split(',')
						if float(x[2]) == 0:		# No FJ - Et, Et_le, Et_ri
							data_Mat[j] = ([float(x[7]),float(x[26]),float(x[27])])
						else:
							data_Mat[j] = ([float(x[11]),float(x[26]),float(x[27])])
						j = j + 1
					c = c + 1
			
			# Number of elements of length = 1 segment
			nELpSeg = nEL/(L/150.)				# Number of elements per segment
			if nELpSeg != round(nELpSeg,0):		# Exit if elementsize is not consistent with segment length
				error_exit = 1
				if L/150. != round(L/150.,0):
					print('Error: L is not a multiple of 150 mm')
				else:
					print('Error: nEL is not suitable')
				break
			else:
				nELpSeg = int(nELpSeg)
			
			# Pick other material values from material table
			if Material_grade == 'user':
				myMaterial = Material_grade
			else:
				myMaterial = Material_grade+'_'+Fraktile
			for intM in range(len(Materialdb)):
				row = Materialdb[intM]
				if row[0] == myMaterial:
					break
					
			# Loop over beam height
			j = 0
			k = 1
			for intH in range(int(H/Hlam)):
				# Loop over beam length
				for intL in range(int(L/150*nELpSeg)):
				# Loop over beam width
					for intW in range(2):
						# Assign material values
						if intW == 0:
							E_1		= data_Mat[j,1]/data_Mat[j,0]*row[1]		#[N/mm^2]
						else:
							E_1		= data_Mat[j,2]/data_Mat[j,0]*row[1]		#[N/mm^2]							
						E_2		= row[2]			#[N/mm^2]
						E_3		= row[3]			#[N/mm^2]
						v_12	= row[4]			#[-]
						v_13	= row[5]			#[-]
						v_23	= row[6]			#[-]			
						G_23	= row[7]			#[N/mm^2]
						G_12	= row[8]			#[N/mm^2]
						G_13	= row[9]			#[N/mm^2]
						f_c0	= row[10]			#[N/mm^2]
						# Timber at supports - compression plasticising, compression and tensile failure are neglected
						if intL < int(nEL_supp_el_0) or intL >= nEL-int(nEL_supp_el_L):
							f_c0= 1000				#[N/mm^2]
						else:
							f_c0= row[10]			#[N/mm^2]	
						# Timber at supports or load introduction - tensile failure is neglected
						if intL < int(nEL_supp_el_0) or intL >= nEL-int(nEL_supp_el_L) or (L_screws_load > 0.0 and Step_n_My != 0 and Moment == 'Qz_3P_area' and intL >= int(nEL/2)-int(nEL_screws_load/2) and intL < int(nEL/2)+int(nEL_screws_load/2) and intH < max(1,int(L_screws_load*nEH/H+1))):
							f_t0= 1000				#[N/mm^2]
						else:
							f_t0= row[11]			#[N/mm^2]
						f_v12	= row[12]			#[N/mm^2
						f_v13	= row[13]			#[N/mm^2]
						f_v23	= row[14]			#[N/mm^2]
						v_21	= (v_12*E_2)/E_1
						v_31	= (v_13*E_3)/E_1
						v_32	= (v_23*E_3)/E_2
						Delta	= 1/(1-v_12*v_21-v_13*v_31-v_23*v_32-2*v_13*v_21*v_32)
						if intW == 1:
							if k == nELpSeg:
								j = j + 1
								k = 1
								if intL == (int(L/150*nELpSeg)-1) and newLam == 1:
									j = j + 54
							else:
								k = k + 1

						# Calculate stiffness matrix for direct input
						D1111 = round((1-v_32*v_23)*Delta*E_1)
						D1122 = round((v_13*v_32+v_12)*Delta*E_2)
						D2222 = round((1-v_13*v_31)*Delta*E_2)
						D1133 = round((v_12*v_23+v_13)*Delta*E_3)
						D2233 = round((v_21*v_13+v_23)*Delta*E_3)
						D3333 = round((1-v_12*v_21)*Delta*E_3)
						D1112 = 0.0
						D2212 = 0.0
						D3312 = 0.0
						D1212 = round(G_12)
						D1113 = 0.0
						D2213 = 0.0
						D3313 = 0.0
						D1213 = 0.0
						D1313 = round(G_13)
						D1123 = 0.0
						D2223 = 0.0
						D3323 = 0.0
						D1223 = 0.0
						D1323 = 0.0
						D2323 = round(G_23)
						
						# Define material models in Abaqus
						# Timber in tension
						myMaterialGL = myModel.Material(name='Timber_H'+str(intH+1)+'_L'+str(intL+1)+'_W'+str(intW+1))
						if Material_law == 'UMAT':
							myMaterialGL.UserMaterial(type=MECHANICAL, mechanicalConstants=(E_1, E_2, E_3, v_12, v_13, v_23, G_12, G_13, G_23, f_c0pl, -f_c0, e_c0pl, E_1pl, f_v12pl1, f_v12pl2*f_v12, e_v21pl, G_12pl, f_v13pl1, f_v13pl2*f_v13, e_v31pl, G_13pl, f_v23pl1, f_v23pl2*f_v23, e_v23pl, G_23pl, f_v12, f_t0))
							if elem_deletion == 1:
								myMaterialGL.Depvar(n=62, deleteVar=9)
							else:
								myMaterialGL.Depvar(n=62)
						else:	
							myMaterialGL.Elastic(type=ANISOTROPIC,table=((D1111, D1122, D2222, D1133, D2233, D3333, D1112, D2212, D3312, D1212, D1113, D2213, D3313, D1213, D1313, D1123, D2223, D3323, D1223, D1323, D2323), ))
							if Plasticising == 'elpl' and Material_law == 'VM':
								myMaterialGL.Plastic(table=((((f_c0pl+(1-f_c0pl)*sin(0.0*pi/2)))*f_c0, e_c0pl*(1-cos(0.0*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.1*pi/2)))*f_c0, e_c0pl*(1-cos(0.1*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.2*pi/2)))*f_c0, e_c0pl*(1-cos(0.2*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.3*pi/2)))*f_c0, e_c0pl*(1-cos(0.3*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.4*pi/2)))*f_c0, e_c0pl*(1-cos(0.4*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.5*pi/2)))*f_c0, e_c0pl*(1-cos(0.5*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.6*pi/2)))*f_c0, e_c0pl*(1-cos(0.6*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.7*pi/2)))*f_c0, e_c0pl*(1-cos(0.7*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.8*pi/2)))*f_c0, e_c0pl*(1-cos(0.8*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(0.9*pi/2)))*f_c0, e_c0pl*(1-cos(0.9*pi/2))*f_c0/E_1), 
															(((f_c0pl+(1-f_c0pl)*sin(1.0*pi/2)))*f_c0, e_c0pl*(1-cos(1.0*pi/2))*f_c0/E_1), 
															(f_c0+0.25*E_1pl, (f_c0/E_1+0.25))))
		
		''' Steel '''
	elif Material == 'Steel':
		myMaterialGL = myModel.Material(name='Steel')
		myMaterialGL.Elastic(table=((210000.0, 0.3), ))
		f_t0 = 235
	
		''' Concrete '''
	elif Material == 'Concrete':
			myMaterialGL = myModel.Material(name='Concrete')
			myMaterialGL.Elastic(table=((33000.0, 0.2), ))	
			f_t0 = 30
	
	''' Material model of the equivalent beam of the bracing '''	 
	if Brace_type == 'springs+beam':
		myMaterialGL = myModel.Material(name='Mat_EquBeamBracing')
		myMaterialGL.Elastic(table=((E_beam_brace, 0.3), ))
		
	if error_exit == 1:
		print('Fatal error')
		break
	################################################################################################################################################################################################################################################
	# the following code does not have to be modified any more
	################################################################################################################################################################################################################################################
	''' Computation of further input parameters ''' 
	################################################################################################################################################################################################################################################
	''' Beam heights ''' 
	########################################
	H_0 = (lowchord[0]-upchord[0])			# Beam height at x = 0
	H_L = (lowchord[-1]-upchord[-1])		# Beam height at x = L
	z_cen_0 = (lowchord[0]+upchord[0])/2	# z-coordinate of the beam centroid at x = 0
	z_cen_L = (lowchord[-1]+upchord[-1])/2	# z-coordinate of the beam centroid at x = L
	
	# z-coordinate of the beam centroid at the support close to x = 0 and L
	nEL_Supp_0 = Pos_supp_0 / sizeEL 
	nEL_Supp_L = Pos_supp_L / sizeEL 
	if nEL_Supp_0 == int(nEL_Supp_0):
		z_Supp_ideal_0 = (lowchord[int(nEL_Supp_0)]+upchord[int(nEL_Supp_0)])/2
		z_Supp_ideal_L = (lowchord[int(nEL_Supp_L)]+upchord[int(nEL_Supp_L)])/2
		z_Supp_real_0 = lowchord[int(nEL_Supp_0)]
		z_Supp_real_L = lowchord[int(nEL_Supp_L)]
	else:	
		z_Supp_ideal_0 = (lowchord[int(nEL_Supp_0)]+upchord[int(nEL_Supp_0)])/2 + ((lowchord[int(nEL_Supp_0+1)]+upchord[int(nEL_Supp_0+1)])/2-(lowchord[int(nEL_Supp_0)]+upchord[int(nEL_Supp_0)])/2) * (nEL_Supp_0-int(nEL_Supp_0))
		z_Supp_ideal_L = (lowchord[int(nEL_Supp_L)]+upchord[int(nEL_Supp_L)])/2 + ((lowchord[int(nEL_Supp_L+1)]+upchord[int(nEL_Supp_L+1)])/2-(lowchord[int(nEL_Supp_L)]+upchord[int(nEL_Supp_L)])/2) * (nEL_Supp_0-int(nEL_Supp_0))
		z_Supp_real_0 = lowchord[int(nEL_Supp_0)] + (lowchord[int(nEL_Supp_0+1)]-lowchord[int(nEL_Supp_0)]) * (nEL_Supp_0-int(nEL_Supp_0))
		z_Supp_real_L = lowchord[int(nEL_Supp_L)] + (lowchord[int(nEL_Supp_L+1)]-lowchord[int(nEL_Supp_L)]) * (nEL_Supp_L-int(nEL_Supp_L))
	
	################################################################################################################################################################################################################################################
	''' Modelling ''' 
	################################################################################################################################################################################################################################################
	''' Parts '''
	########################################
	''' Nodes names '''
	# A list of the node names is generated
	myNodeLabels = []
	[ ( myNodeLabels.append(int(10010010 + intH*10 + intB*10000 + intL*10000000))) for intL in range(nEL+1) for intH in range(nEH+1) for intB in range(nEB+1) ]	# Corner nodes of C3D20R elements
	[ ( myNodeLabels.append(int(15010010 + intH*10 + intB*10000 + intL*10000000))) for intL in range(nEL  ) for intH in range(nEH+1) for intB in range(nEB+1) ]	# Intermediate nodes over the length of the C3D20R elements
	[ ( myNodeLabels.append(int(10015010 + intH*10 + intB*10000 + intL*10000000))) for intL in range(nEL+1) for intH in range(nEH+1) for intB in range(nEB  ) ]	# Intermediate nodes over the width of the C3D20R elements
	[ ( myNodeLabels.append(int(10010015 + intH*10 + intB*10000 + intL*10000000))) for intL in range(nEL+1) for intH in range(nEH  ) for intB in range(nEB+1) ]	# Intermediate nodes over the height of the C3D20R elements
	if Brace_type == 'springs+beam':	# Nodes of the equivalent beam of the bracing
		myNodeLabelsBracing = []
		[ ( myNodeLabelsBracing.append(int(intL+1))) for intL in range(nEL+1) ]

	''' Nodes coordinates '''
	# A list of the node coordinates is generated
	Coordinates = []
	[Coordinates.append((sizeEL*intL, 		sizeEB*intB, 				upchord[intL] + intH * (lowchord[intL] - upchord[intL])/ nEH)) 																					for intL in range(nEL+1) 	for intH in range(nEH+1) 	for intB in range(-int(nEB/2), int(nEB/2)+1) 	]
	[Coordinates.append((sizeEL*(intL+0.5), sizeEB*intB, 				(upchord[intL] + upchord[intL + 1]) / 2 + intH * ((lowchord[intL] + lowchord[intL + 1]) / 2 - (upchord[intL] + upchord[intL + 1]) / 2)/ nEH)) 	for intL in range(nEL) 		for intH in range(nEH+1) 	for intB in range(-int(nEB/2), int(nEB/2)+1) 	]
	[Coordinates.append((sizeEL*intL, 		sizeEB*0.5+intB*sizeEB, 	upchord[intL] + intH * (lowchord[intL] - upchord[intL])/ nEH)) 																					for intL in range(nEL+1) 	for intH in range(nEH+1) 	for intB in range(-int(nEB/2), int(nEB/2))	]
	[Coordinates.append((sizeEL*intL, 		sizeEB*intB, 				upchord[intL] + (intH + 0.5) * (lowchord[intL] - upchord[intL])/ nEH)) 																			for intL in range(nEL+1) 	for intH in range(nEH) 		for intB in range(-int(nEB/2), int(nEB/2)+1) 	]
	if Brace_type == 'springs+beam':	# Nodes of the equivalent beam of the bracing
		myNodesCoordsBracing = []
		[myNodesCoordsBracing.append((sizeEL*intL, 	1000, 	(t_up1+t_u1)/2 + ((t_up2+t_u2)/2-(t_up1+t_u1)/2)/L * intL*sizeEL)) 	for intL in range(nEL+1) ]
	
		''' Coordinates of the imperfections at each node'''	
	# LBA - no imperfections
	if Cal_type == 'LBA':	
		L1 = np.array(Coordinates)
		myNodesCoords = L1.tolist()
	
	# GMNIA - imperfections
	elif Cal_type == 'GMNIA':	
		# Imperfections
		if Imp_type == 'equi' or Imp_type == 'meas':
			# bow
			bow = []
			[bow.append((0.0, ey[intL], ez[intL])) 		for intL in range(nEL+1) 	for intH in range(nEH+1) 	for intB in range(-int(nEB/2), int(nEB/2)+1)	]		# Corner nodes
			[bow.append((0.0, ey_IN[intL], ez_IN[intL])) 	for intL in range(nEL) 		for intH in range(nEH+1) 	for intB in range(-int(nEB/2), int(nEB/2)+1)	]		# Intermediate nodes lenght
			[bow.append((0.0, ey[intL], ez[intL])) 		for intL in range(nEL+1) 	for intH in range(nEH+1) 	for intB in range(-int(nEB/2), int(nEB/2))	]		# Intermediate nodes width
			[bow.append((0.0, ey[intL], ez[intL])) 		for intL in range(nEL+1) 	for intH in range(nEH) 		for intB in range(-int(nEB/2), int(nEB/2)+1)	]		# Intermediate nodes height
					
			# Correction of the z coordinates, so that z = 0 is in thr cross-section centroid. Otherwise the rotation forumula would be a lot more complicated.
			z_correction = np.add(lowchord, -np.add(-np.array(upchord), lowchord)/2)
			upchord_Czero = np.add(upchord, -z_correction)
			lowchord_Czero = np.add(lowchord, -z_correction)
			
			# twist
			twist = []	
			[twist.append((sizeEL*intL, 		sizeEB * intB * math.cos(etheta[intL]) + (intH*(lowchord_Czero[intL]-upchord_Czero[intL])/nEH+upchord_Czero[intL]) * -math.sin(etheta[intL]) , sizeEB * intB * math.sin(etheta[intL]) + (intH*(lowchord_Czero[intL]-upchord_Czero[intL])/nEH+upchord_Czero[intL]) * math.cos(etheta[intL]) + z_correction[intL]))																																																			for intL in range(nEL+1) 	for intH in range(nEH+1) 	for intB in range(-int(nEB/2), int(nEB/2)+1)	]
			[twist.append((sizeEL*(intL+0.5), 	sizeEB * intB * math.cos(etheta_IN[intL]) + (intH*((lowchord_Czero[intL]+lowchord_Czero[intL+1])/2-(upchord_Czero[intL]+upchord_Czero[intL+1])/2)/nEH+(upchord_Czero[intL]+upchord_Czero[intL+1])/2)	* -math.sin(etheta_IN[intL]) , sizeEB * intB * math.sin(etheta_IN[intL]) + (intH*((lowchord_Czero[intL]+lowchord_Czero[intL+1])/2-(upchord_Czero[intL]+upchord_Czero[intL+1])/2)/nEH+(upchord_Czero[intL]+upchord_Czero[intL+1])/2) * math.cos(etheta_IN[intL]) + (z_correction[intL]+z_correction[intL+1])/2 )) 	for intL in range(nEL) 		for intH in range(nEH+1) 	for intB in range(-int(nEB/2), int(nEB/2)+1)	]
			[twist.append((sizeEL*intL	, 		sizeEB * (intB+0.5) * math.cos(etheta[intL]) + (intH*(lowchord_Czero[intL]-upchord_Czero[intL])/nEH+upchord_Czero[intL]) * -math.sin(etheta[intL]) , sizeEB * (intB+0.5) * math.sin(etheta[intL]) + (intH*(lowchord_Czero[intL]-upchord_Czero[intL])/nEH+upchord_Czero[intL]) * math.cos(etheta[intL]) + z_correction[intL]))																																																for intL in range(nEL+1) 	for intH in range(nEH+1) 	for intB in range(-int(nEB/2), int(nEB/2))	]
			[twist.append((sizeEL*intL, 		sizeEB * intB * math.cos(etheta[intL]) + ((intH+0.5)*(lowchord_Czero[intL]-upchord_Czero[intL])/nEH+upchord_Czero[intL]) * -math.sin(etheta[intL]) , sizeEB * intB * math.sin(etheta[intL]) + ((intH+0.5)*(lowchord_Czero[intL]-upchord_Czero[intL])/nEH+upchord_Czero[intL]) * math.cos(etheta[intL]) + z_correction[intL]))																																																for intL in range(nEL+1) 	for intH in range(nEH) 		for intB in range(-int(nEB/2), int(nEB/2)+1)	]

			# Definition of lists of node coordinates and imperfection values
			L1 = np.array(Coordinates)
			L2 = np.array(bow)
			L3 = np.array(twist)

			# Add imperfection to ideal beam coordinates for getting the final beam geometry
			if Imp_shape == 'bow':
				L4 = np.add(L1, L2)
				myNodesCoords = L4.tolist()
			elif Imp_shape == 'twist':
				myNodesCoords = L3.tolist()
			elif Imp_shape == 'bow+twist':
				L5 = np.add(L2, L3)
				myNodesCoords = L5.tolist()
				
		# No imperfections
		else:
			L1 = np.array(Coordinates)
			myNodesCoords = L1.tolist()
	
	''' Elements '''
	# A list of the element names is generated
	myElementsLabel = []
	[ ( myElementsLabel.append(int(10010010 + intH*10 + intB*10000 + intL*10000000))) for intL in range(nEL) for intH in range(nEH) for intB in range(nEB) ]
	if Brace_type == 'springs+beam':	# Equivalent beam bracing
		myElementsLabelBracing = []
		[ ( myElementsLabelBracing.append(int(intL+1))) for intL in range(nEL)]
	
	# A list of the node names that are part of each element is generated	
	myElementConns = []
	[myElementConns.append((int(10010010 + intB*10000 + intH*10 + intL*10000000), int(10020010 + intB*10000 + intH*10 + intL*10000000), int(10020020 + intB*10000 + intH*10 + intL*10000000), int(10010020 + intB*10000 + intH*10 + intL*10000000),int(20010010 + intB*10000 + intH*10 + intL*10000000), int(20020010 + intB*10000 + intH*10 + intL*10000000), int(20020020 + intB*10000 + intH*10 + intL*10000000), int(20010020 + intB*10000 + intH*10 + intL*10000000), int(10015010 + intB*10000 + intH*10 + intL*10000000), int(10020015 + intB*10000 + intH*10 + intL*10000000), int(10015020 + intB*10000 + intH*10 + intL*10000000), int(10010015 + intB*10000 + intH*10 + intL*10000000), int(20015010 + intB*10000 + intH*10 + intL*10000000), int(20020015 + intB*10000 + intH*10 + intL*10000000), int(20015020 + intB*10000 + intH*10 + intL*10000000), int(20010015 + intB*10000 + intH*10 + intL*10000000), int(15010010 + intB*10000 + intH*10 + intL*10000000), int(15020010 + intB*10000 + intH*10 + intL*10000000), int(15020020 + intB*10000 + intH*10 + intL*10000000), int(15010020 + intB*10000 + intH*10 + intL*10000000),)) for intL in range(nEL) for intH in range(nEH) for intB in range(nEB) ]
	if Brace_type == 'springs+beam':	# Equivalent beam bracing
		myElementConnsBracing = []
		[myElementConnsBracing.append((int(intL+1), int(intL+2))) for intL in range(nEL)]
		
	# Create Part
	myPartbeam = myModel.PartFromNodesAndElements(name='beam', dimensionality=THREE_D, type=DEFORMABLE_BODY, nodes=(myNodeLabels, myNodesCoords), elements=(("C3D20R", myElementsLabel, myElementConns),))
	if Brace_type == 'springs+beam':	# Equivalent beam bracing
		myPartEquBeamBracing = myModel.PartFromNodesAndElements(name='EquBeamBracing', dimensionality=THREE_D, type=DEFORMABLE_BODY, nodes=(myNodeLabelsBracing, myNodesCoordsBracing), elements=(("B31", myElementsLabelBracing, myElementConnsBracing),))
	
	########################################
	''' Sections '''
	########################################
	if KaReMo == 0 and KaReMo_y == 0:
		# Tension zone
		myElements = []
		[myElements.append(int(10010010+intL*10000000+intB*10000+intH*10), ) for intL in range(nEL) for intB in range(nEB) for intH in range(int(nEH/2),nEH-nEdgelam_Tens)]
		TensZone = myPartbeam.elements.sequenceFromLabels(labels=myElements)
		myPartbeam.Set(elements=TensZone, name='TensZone')
		TensZone = myPartbeam.sets['TensZone']
		if Material == 'Timber':
			mySectTensZone = myModel.HomogeneousSolidSection(name='TensZone', material='Timber_tens')
		elif Material == 'Steel':
			mySectTensZone = myModel.HomogeneousSolidSection(name='TensZone', material='Steel')
		elif Material == 'Concrete':
			mySectTensZone = myModel.HomogeneousSolidSection(name='TensZone', material='Concrete')	
		myPartbeam.SectionAssignment(region=TensZone, sectionName='TensZone')
		
		# Compression zone
		myElements = []
		[myElements.append(int(10010010+intL*10000000+intB*10000+intH*10), ) for intL in range(nEL) for intB in range(nEB) for intH in range(nEdgelam_Comp,int(nEH/2))]
		CompZone = myPartbeam.elements.sequenceFromLabels(labels=myElements)

		myPartbeam.Set(elements=CompZone, name='CompZone')
		CompZone = myPartbeam.sets['CompZone']
		if Material == 'Timber':
			mySectCompZone = myModel.HomogeneousSolidSection(name='CompZone', material='Timber_comp')
		elif Material == 'Steel':
			mySectCompZone = myModel.HomogeneousSolidSection(name='CompZone', material='Steel')
		elif Material == 'Concrete':
			mySectCompZone = myModel.HomogeneousSolidSection(name='CompZone', material='Concrete')	
		myPartbeam.SectionAssignment(region=CompZone, sectionName='CompZone')
		
		# Load introduction where tensile failure is neglected
		if L_screws_load > 0.0 and nEL_screws_load != 0 and Step_n_My != 0 and (Moment == 'Qz_3P_area' or Moment == 'Exp_LTB'):
			if (nEL_screws_load % 2) == 0:		# check if even
				minEL = int(round(Pos_Qz_3P/sizeEL,0)-nEL_screws_load/2)
				maxEL = int(round(Pos_Qz_3P/sizeEL,0)+nEL_screws_load/2)
			else:
				if Pos_Qz_3P > sizeEL*int(round(Pos_Qz_3P/sizeEL,0)):
					minEL = int(round(Pos_Qz_3P/sizeEL,0) - (nEL_screws_load-1)/2)
					maxEL = int(round(Pos_Qz_3P/sizeEL,0) + 1 + (nEL_screws_load-1)/2)
				else:
					minEL = int(round(Pos_Qz_3P/sizeEL,0) - 1 - (nEL_screws_load-1)/2)
					maxEL = int(round(Pos_Qz_3P/sizeEL,0) + (nEL_screws_load-1)/2)
		
			myElements = []
			[myElements.append(int(10010010+intL*10000000+intB*10000+intH*10), ) for intL in range(minEL,maxEL) for intB in range(nEB) for intH in range(max(1,int(L_screws_load*nEH/H+1)))]
			LoadZone = myPartbeam.elements.sequenceFromLabels(labels=myElements)
			myPartbeam.Set(elements=LoadZone, name='LoadZone')
			LoadZone = myPartbeam.sets['LoadZone']
			if Material == 'Timber':
				mySectLoadZone = myModel.HomogeneousSolidSection(name='LoadZone', material='Timber_load')
			elif Material == 'Steel':
				mySectLoadZone = myModel.HomogeneousSolidSection(name='LoadZone', material='Steel')
			elif Material == 'Concrete':
				mySectLoadZone = myModel.HomogeneousSolidSection(name='LoadZone', material='Concrete')	
			myPartbeam.SectionAssignment(region=LoadZone, sectionName='LoadZone')

		# Support zone where compression plasticising and compression and tensile failure are neglected
		if nEL_supp_el_0 != 0 or nEL_supp_el_L != 0:
			myElements = []
			[myElements.append(int(10010010+intL*10000000+intB*10000+intH*10), ) for intL in range(int(nEL_supp_el_0)) for intB in range(nEB) for intH in range(nEH)]
			[myElements.append(int(10010010+intL*10000000+intB*10000+intH*10), ) for intL in range((nEL-int(nEL_supp_el_L)),nEL) for intB in range(nEB) for intH in range(nEH)]
			SuppZone = myPartbeam.elements.sequenceFromLabels(labels=myElements)
			myPartbeam.Set(elements=SuppZone, name='SuppZone')
			SuppZone = myPartbeam.sets['SuppZone']
			if Material == 'Timber':
				mySectSuppZone = myModel.HomogeneousSolidSection(name='SuppZone', material='Timber_supp')
			elif Material == 'Steel':
				mySectSuppZone = myModel.HomogeneousSolidSection(name='SuppZone', material='Steel')
			elif Material == 'Concrete':
				mySectSuppZone = myModel.HomogeneousSolidSection(name='SuppZone', material='Concrete')	
			myPartbeam.SectionAssignment(region=SuppZone, sectionName='SuppZone')

		# Edge lamellae tension
		if nEdgelam_Tens != 0:
			myElements = []
			[myElements.append(int(10010010+intL*10000000+intB*10000+intH*10), ) for intL in range(nEL) for intB in range(nEB) for intH in range(nEH-nEdgelam_Tens,nEH)]
			Edgelam_Tens = myPartbeam.elements.sequenceFromLabels(labels=myElements)
			myPartbeam.Set(elements=Edgelam_Tens, name='Edgelam_Tens')
			Edgelam_Tens = myPartbeam.sets['Edgelam_Tens']
			mySectEdgelam_Tens = myModel.HomogeneousSolidSection(name='Edgelam_Tens', material='Timber_Edgelam_Tens')
			myPartbeam.SectionAssignment(region=Edgelam_Tens, sectionName='Edgelam_Tens')
		
		# Edge lamellae compression
		if nEdgelam_Comp != 0:
			myElements = []
			[myElements.append(int(10010010+intL*10000000+intB*10000+intH*10), ) for intL in range(nEL) for intB in range(nEB) for intH in range(nEdgelam_Comp)]
			Edgelam_Comp = myPartbeam.elements.sequenceFromLabels(labels=myElements)
			myPartbeam.Set(elements=Edgelam_Comp, name='Edgelam_Comp')
			Edgelam_Comp = myPartbeam.sets['Edgelam_Comp']
			mySectEdgelam_Comp = myModel.HomogeneousSolidSection(name='Edgelam_Comp', material='Timber_Edgelam_Comp')
			myPartbeam.SectionAssignment(region=Edgelam_Comp, sectionName='Edgelam_Comp')
			
	else:
		# Loop over the beam height
		for intH in range(int(H/Hlam)):
			# Loop over beam length
			for intL in range(int(L/150*nELpSeg)):
				# Loop over beam width
				for intW in range(2):
					myElements = []
					if intW == 0:
						[myElements.append(int(10010010+intL*10000000+intB*10000+intH*10), ) for intB in range(int(nEB/2))]
					else:
						[myElements.append(int(10010010+intL*10000000+intB*10000+intH*10), ) for intB in range(int(nEB/2),nEB)]
					myElem = myPartbeam.elements.sequenceFromLabels(labels=myElements)
					myName='Timber_H'+str(intH+1)+'_L'+str(intL+1)+'_W'+str(intW+1)
					myPartbeam.Set(elements=myElem, name=myName)
					mySet = myPartbeam.sets[myName]
					mySect = myModel.HomogeneousSolidSection(name=myName, material=myName)
					myPartbeam.SectionAssignment(region=mySet, sectionName=myName)
		
	# Equivalent beam bracing
	if Brace_type == 'springs+beam':
		# Create cross-section
		myModel.RectangularProfile(name='CS_EquBeamBracing', a = 100, b = B_beam_brace)		# a = in z direction, b = in y direction
		myModel.BeamSection(name='EquBeamBracing', integration=DURING_ANALYSIS, poissonRatio=0.0, profile='CS_EquBeamBracing', material='Mat_EquBeamBracing', temperatureVar=LINEAR, consistentMassMatrix=False)
		EquBeamBracing = myPartEquBeamBracing.elements.sequenceFromLabels(labels=myElementsLabelBracing)
		myPartEquBeamBracing.Set(elements=EquBeamBracing, name='EquBeamBracing')
		EquBeamBracing = myPartEquBeamBracing.sets['EquBeamBracing']		
		myPartEquBeamBracing.SectionAssignment(region=EquBeamBracing, sectionName='EquBeamBracing')
		myModel.parts['EquBeamBracing'].assignBeamSectionOrientation(region=EquBeamBracing, method=N1_COSINES, n1=(0.0, 0.0, -1.0))	# Assign beam orientation
		
	########################################
	''' Material Orientation '''
	########################################
	# The material orientation is chosen depending on the lower edge of the beam, perpendicular lower edge
	
	# Chose surfaces on the lower edge
	myElements = []
	[myElements.append(int(10010010+intL*10000000+intB*10000+(nEH-1)*10), ) for intL in range(nEL) for intB in range(nEB) ]
	elemLowEdge = myPartbeam.elements.sequenceFromLabels(labels=myElements)
	myPartbeam.Set(elements=elemLowEdge, name='Elements_LowEdge')	
	myPartbeam.Surface(face5Elements=elemLowEdge, name='LowEdge')
	
	if Material == 'Timber':
		# Material orientation depending on the lower edge
		region = regionToolset.Region(elements = mdb.models[modelName].parts['beam'].elements)	
		normalAxisRegion = mdb.models[modelName].parts['beam'].surfaces['LowEdge']
		if Mat_axis_vert == 'R':	# long = L, vert = R
			mdb.models[modelName].parts['beam'].MaterialOrientation(
				region=region, orientationType=DISCRETE, axis=AXIS_1, 
				normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
				flipNormalDirection=False, normalAxisDirection=AXIS_3, 
				primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 1.0, 0.0),
				primaryAxisDirection=AXIS_2, flipPrimaryDirection=False, 
				additionalRotationType=ROTATION_NONE, angle=90, 
				additionalRotationField='', stackDirection=STACK_3)
		else:						# long = L, vert = T
			mdb.models[modelName].parts['beam'].MaterialOrientation(
				region=region, orientationType=DISCRETE, axis=AXIS_1, 
				normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
				flipNormalDirection=False, normalAxisDirection=AXIS_3, 
				primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 1.0, 0.0),
				primaryAxisDirection=AXIS_2, flipPrimaryDirection=False, 
				additionalRotationType=ROTATION_NONE, angle=0, 
				additionalRotationField='', stackDirection=STACK_3)
					
	#######################################
	''' Assembly '''
	########################################
	myAssembly = myModel.rootAssembly
	# Beam
	myInstanceBeam = myAssembly.Instance(name='instbeam', part=myPartbeam, dependent=ON)
	# Equivalent beam bracing
	if Brace_type == 'springs+beam':
		myInstanceEquBeamBracing = myAssembly.Instance(name='instEquBeamBracing', part=myPartEquBeamBracing, dependent=ON)
		
	########################################
	''' Reference nodes '''
	########################################
	''' Supports '''
	# Ideal supports at cross-section centroid
	if Supp_type == 'ideal':
		# Apply load Nx eccentricity by moving supports
		if Step_n_Nx != 0 and Load_type_Nx == 'w':
			# Support at x = 0
			myAssembly.ReferencePoint(point=(0.0, ey_Nx_0, z_Supp_ideal_0+ez_Nx_0))
			myAssembly.features.changeKey(fromName='RP-1', toName='RP_Supp_0')
			RP_Supp_0 = (myAssembly.referencePoints.findAt((0.0, ey_Nx_0, z_Supp_ideal_0+ez_Nx_0)),)
			RP_Supp_0 = myAssembly.Set(referencePoints=RP_Supp_0, name='RP_Supp_0')
			# Support at x = L
			myAssembly.ReferencePoint(point=(L, ey_Nx_L, z_Supp_ideal_L+ez_Nx_L))
			myAssembly.features.changeKey(fromName='RP-1', toName='RP_Supp_L')
			RP_Supp_L = (myAssembly.referencePoints.findAt((L, ey_Nx_L, z_Supp_ideal_L+ez_Nx_L)),)
			RP_Supp_L = myAssembly.Set(referencePoints=RP_Supp_L, name='RP_Supp_L')
		else:
			# Support at x = 0
			myAssembly.ReferencePoint(point=(0.0, 0.0, z_Supp_ideal_0))
			myAssembly.features.changeKey(fromName='RP-1', toName='RP_Supp_0')
			RP_Supp_0 = (myAssembly.referencePoints.findAt((0.0, 0.0, z_Supp_ideal_0)),)
			RP_Supp_0 = myAssembly.Set(referencePoints=RP_Supp_0, name='RP_Supp_0')
			# Support at x = L
			myAssembly.ReferencePoint(point=(L, 0.0, z_Supp_ideal_L))
			myAssembly.features.changeKey(fromName='RP-1', toName='RP_Supp_L')
			RP_Supp_L = (myAssembly.referencePoints.findAt((L, 0.0, z_Supp_ideal_L)),)
			RP_Supp_L = myAssembly.Set(referencePoints=RP_Supp_L, name='RP_Supp_L')
		
	# Realistic supporting conditions at lower edge of the beam
	elif Supp_type == 'real':
		# Support at x = 0
		myAssembly.ReferencePoint(point=(Pos_supp_0, 0.0, z_Supp_real_0))
		myAssembly.features.changeKey(fromName='RP-1', toName='RP_Supp_0')
		RP_Supp_0 = (myAssembly.referencePoints.findAt((Pos_supp_0, 0.0, z_Supp_real_0)),)
		RP_Supp_0 = myAssembly.Set(referencePoints=RP_Supp_0, name='RP_Supp_0')
		
		# Support at x = L
		myAssembly.ReferencePoint(point=(Pos_supp_L, 0.0, z_Supp_real_L))
		myAssembly.features.changeKey(fromName='RP-1', toName='RP_Supp_L')
		RP_Supp_L = (myAssembly.referencePoints.findAt((Pos_supp_L, 0.0, z_Supp_real_L)),)
		RP_Supp_L = myAssembly.Set(referencePoints=RP_Supp_L, name='RP_Supp_L')

		# Fork support at x = 0
		myAssembly.ReferencePoint(point=(Pos_supp_0, 0.0, z_Supp_ideal_0))
		myAssembly.features.changeKey(fromName='RP-1', toName='RP_ForkSupp_0')
		RP_ForkSupp_0 = (myAssembly.referencePoints.findAt((Pos_supp_0, 0.0, z_Supp_ideal_0)),)
		RP_ForkSupp_0 = myAssembly.Set(referencePoints=RP_ForkSupp_0, name='RP_ForkSupp_0')

		# Fork support at x = L
		myAssembly.ReferencePoint(point=(Pos_supp_L, 0.0, z_Supp_ideal_L))
		myAssembly.features.changeKey(fromName='RP-1', toName='RP_ForkSupp_L')
		RP_ForkSupp_L = (myAssembly.referencePoints.findAt((Pos_supp_L, 0.0, z_Supp_ideal_L)),)
		RP_ForkSupp_L = myAssembly.Set(referencePoints=RP_ForkSupp_L, name='RP_ForkSupp_L')
	
	''' Bracing '''	
	if Brace_type == 'springs+beam' or Brace_type == 'NodestoGround':	
		for intA in range(len(pos_brace)):
			mycoordinates = myModel.parts['beam'].nodes[int((round(nEL*pos_brace[intA],0))*(nEH+1)*(nEB+1)+((min_nNH_brace+max_nNH_brace)/2-1)*(nEB+1)+(nEB+1)*0.5)].coordinates
			myAssembly.ReferencePoint(point=mycoordinates)
			myNameRP = 'RP_Brace_'+str(intA)
			myName = 'RP_Brace_'+str(intA)
			myAssembly.features.changeKey(fromName='RP-1', toName=myNameRP)
			if intA == 0:
				RP_Brace_0 = (myAssembly.referencePoints.findAt(mycoordinates),)
				RP_Brace_0 = myAssembly.Set(referencePoints=RP_Brace_0, name=myName)
			elif intA == 1:
				RP_Brace_1 = (myAssembly.referencePoints.findAt(mycoordinates),)
				RP_Brace_1 = myAssembly.Set(referencePoints=RP_Brace_1, name=myName)
			elif intA == 2:
				RP_Brace_2 = (myAssembly.referencePoints.findAt(mycoordinates),)
				RP_Brace_2 = myAssembly.Set(referencePoints=RP_Brace_2, name=myName)
			elif intA == 3:
				RP_Brace_3 = (myAssembly.referencePoints.findAt(mycoordinates),)
				RP_Brace_3 = myAssembly.Set(referencePoints=RP_Brace_3, name=myName)	
			elif intA == 4:
				RP_Brace_4 = (myAssembly.referencePoints.findAt(mycoordinates),)
				RP_Brace_4 = myAssembly.Set(referencePoints=RP_Brace_4, name=myName)	
			elif intA == 5:
				RP_Brace_5 = (myAssembly.referencePoints.findAt(mycoordinates),)
				RP_Brace_5 = myAssembly.Set(referencePoints=RP_Brace_5, name=myName)	
			elif intA == 6:
				RP_Brace_6 = (myAssembly.referencePoints.findAt(mycoordinates),)
				RP_Brace_6 = myAssembly.Set(referencePoints=RP_Brace_6, name=myName)	
			elif intA == 7:
				RP_Brace_7 = (myAssembly.referencePoints.findAt(mycoordinates),)
				RP_Brace_7 = myAssembly.Set(referencePoints=RP_Brace_7, name=myName)	
			elif intA == 8:
				RP_Brace_8 = (myAssembly.referencePoints.findAt(mycoordinates),)
				RP_Brace_8 = myAssembly.Set(referencePoints=RP_Brace_8, name=myName)	
			elif intA == 9:
				RP_Brace_9 = (myAssembly.referencePoints.findAt(mycoordinates),)
				RP_Brace_9 = myAssembly.Set(referencePoints=RP_Brace_9, name=myName)	
	
	''' Actions '''
	# Axial actions
	if Step_n_Nx != 0 and Load_type_Nx == 'w' and Supp_type == 'real':
		# Loading at x = 0
		myAssembly.ReferencePoint(point=(0.0, ey_Nx_0, z_Supp_ideal_0+ez_Nx_0))
		myAssembly.features.changeKey(fromName='RP-1', toName='RP_Nx_0')
		RP_Nx_0 = (myAssembly.referencePoints.findAt((0.0, ey_Nx_0, z_Supp_ideal_0+ez_Nx_0)),)
		RP_Nx_0 = myAssembly.Set(referencePoints=RP_Nx_0, name='RP_Nx_0')
		# Loading at x = L
		myAssembly.ReferencePoint(point=(L, ey_Nx_L, z_Supp_ideal_L+ez_Nx_L))
		myAssembly.features.changeKey(fromName='RP-1', toName='RP_Nx_L')
		RP_Nx_L = (myAssembly.referencePoints.findAt((L, ey_Nx_L, z_Supp_ideal_L+ez_Nx_L)),)
		RP_Nx_L = myAssembly.Set(referencePoints=RP_Nx_L, name='RP_Nx_L')
	
	if Step_n_Nx != 0 and Load_type_Nx == 'Exp_LTB' and Supp_type == 'real':		
		# Loading at x = 0
		myAssembly.ReferencePoint(point=(0.0, ey_Nx_0, z_Supp_ideal_0+ez_Nx_0))
		myAssembly.features.changeKey(fromName='RP-1', toName='RP_Nx_0')
		RP_Nx_0 = (myAssembly.referencePoints.findAt((0.0, ey_Nx_0, z_Supp_ideal_0+ez_Nx_0)),)
		RP_Nx_0 = myAssembly.Set(referencePoints=RP_Nx_0, name='RP_Nx_0')
		# Loading at x = L
		myAssembly.ReferencePoint(point=(L, ey_Nx_L, z_Supp_ideal_L+ez_Nx_L))
		myAssembly.features.changeKey(fromName='RP-1', toName='RP_Nx_L')
		RP_Nx_L = (myAssembly.referencePoints.findAt((L, ey_Nx_L, z_Supp_ideal_L+ez_Nx_L)),)
		RP_Nx_L = myAssembly.Set(referencePoints=RP_Nx_L, name='RP_Nx_L')
		
	# Torsional actions
	if Step_n_Mx != 0:
		# Node numbers at x location of Mx
		nN_Mx = nNH*nNB*int(round(Pos_Mx/sizeEL,0))+int(nEH/2)*nNB+int(nEB/2)+1
			
		# Coordinates of the reference point including effects of imperfections
		Coords_Mx_x = myNodesCoords[nN_Mx-1][0]
		Coords_Mx_y = myNodesCoords[nN_Mx-1][1]
		Coords_Mx_z = myNodesCoords[nN_Mx-1][2] 
			
		myAssembly.ReferencePoint(point=(Coords_Mx_x, Coords_Mx_y, Coords_Mx_z))
		myAssembly.features.changeKey(fromName='RP-1', toName='RP_Mx')
		RP_Mx = (myAssembly.referencePoints.findAt((Coords_Mx_x, Coords_Mx_y, Coords_Mx_z)),)
		RP_Mx = myAssembly.Set(referencePoints=RP_Mx, name='RP_Mx')
	
	# Actions causing My
	if Step_n_My != 0:

		# Point load - 3-point bending
		if Moment == 'Qz_3P': 
			# Node numbers at x location of F
			nN_F_upchord 	= nNH*nNB*int(round(Pos_Qz_3P/sizeEL,0))+int(nEB/2)+1				    # Node number at x location of F in the mid of the upper chord
			nN_F_mid 		= nNH*nNB*int(round(Pos_Qz_3P/sizeEL,0))+int(nEH/2)*nNB+int(nEB/2)+1	# Node number at x location of F in the centroid
			nN_F_lowchord	= nNH*nNB*int(round(Pos_Qz_3P/sizeEL,0))+nEH*nNB+int(nEB/2)+1		    # Node number at x location of F in the mid of the lower chord	
			rot_F = (myNodesCoords[nN_F_upchord-1][1] - myNodesCoords[nN_F_lowchord-1][1])/(myNodesCoords[nN_F_lowchord-1][2] - myNodesCoords[nN_F_upchord-1][2])	# Cross-section rotation at F due to imperfections
			
			# Coordinates of the reference point including effects of imperfections
			Coords_3P_x = myNodesCoords[nN_F_mid-1][0]
			Coords_3P_y = myNodesCoords[nN_F_mid-1][1] + ey_Qz * math.cos(rot_F) - ez_Qz * math.sin(rot_F)
			Coords_3P_z = myNodesCoords[nN_F_mid-1][2] + ez_Qz * math.cos(rot_F) + ey_Qz * math.sin(rot_F) 
			
			# Reference point
			myAssembly.ReferencePoint(point=(Coords_3P_x, Coords_3P_y, Coords_3P_z))
			myAssembly.features.changeKey(fromName='RP-1', toName='RP_Qz_3P')
			RP_Qz_3P = (myAssembly.referencePoints.findAt((Coords_3P_x, Coords_3P_y, Coords_3P_z)),)
			RP_Qz_3P = myAssembly.Set(referencePoints=RP_Qz_3P, name='RP_Qz_3P')
		
		# Point load - 3-point bending - distributed over an area
		if Moment == 'Qz_3P_area' or Moment == 'Exp_LTB': 
			# Node numbers at x location of F
			nN_F_upchord 	= nNH*nNB*int(round(Pos_Qz_3P/sizeEL,0))+int(nEB/2)+1				    # Node number at x location of F in the mid of the upper chord
			nN_F_mid 		= nNH*nNB*int(round(Pos_Qz_3P/sizeEL,0))+int(nEH/2)*nNB+int(nEB/2)+1	# Node number at x location of F in the centroid
			nN_F_lowchord	= nNH*nNB*int(round(Pos_Qz_3P/sizeEL,0))+nEH*nNB+int(nEB/2)+1		    # Node number at x location of F in the mid of the lower chord	
			rot_F = (myNodesCoords[nN_F_upchord-1][1] - myNodesCoords[nN_F_lowchord-1][1])/(myNodesCoords[nN_F_lowchord-1][2] - myNodesCoords[nN_F_upchord-1][2])	# Cross-section rotation at F due to imperfections
			
			Coords_3P_x = Pos_Qz_3P
			Coords_3P_y = myNodesCoords[nN_F_mid-1][1] + ey_Qz * math.cos(rot_F) - ez_Qz * math.sin(rot_F)
			Coords_3P_z = myNodesCoords[nN_F_mid-1][2] + ez_Qz * math.cos(rot_F) + ey_Qz * math.sin(rot_F)	
			
			myAssembly.ReferencePoint(point=(Coords_3P_x, Coords_3P_y, Coords_3P_z))
			myAssembly.features.changeKey(fromName='RP-1', toName='RP_Qz_3P_area')
			RP_Qz_3P_area = (myAssembly.referencePoints.findAt((Coords_3P_x, Coords_3P_y, Coords_3P_z)),)
			RP_Qz_3P_area = myAssembly.Set(referencePoints=RP_Qz_3P_area, name='RP_Qz_3P_area')
		
		# 2 point loads - 4-point bending
		if Moment == 'Qz_4P' or Moment == 'Qz_4P_area_T': 
			# Load closer to x = 0
			# Node numbers at x location of F
			nN_F_upchord 	= nNH*nNB*int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0))+int(nEB/2)+1				    # Node number at x location of F in the mid of the upper chord
			nN_F_mid 		= nNH*nNB*int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0))+int(nEH/2)*nNB+int(nEB/2)+1	# Node number at x location of F in the centroid
			nN_F_lowchord	= nNH*nNB*int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0))+nEH*nNB+int(nEB/2)+1		    # Node number at x location of F in the mid of the lower chord	
			rot_F = (myNodesCoords[nN_F_upchord-1][1] - myNodesCoords[nN_F_lowchord-1][1])/(myNodesCoords[nN_F_lowchord-1][2] - myNodesCoords[nN_F_upchord-1][2])	# Cross-section rotation at F due to imperfections
			
			# Coordinates of the reference point including effects of imperfections
			Coords_4P_0_x = myNodesCoords[nN_F_mid-1][0]
			Coords_4P_0_y = myNodesCoords[nN_F_mid-1][1] + ey_Qz * math.cos(rot_F) - ez_Qz * math.sin(rot_F)
			Coords_4P_0_z = myNodesCoords[nN_F_mid-1][2] + ez_Qz * math.cos(rot_F) + ey_Qz * math.sin(rot_F) 
			
			# Reference point
			myAssembly.ReferencePoint(point=(Coords_4P_0_x, Coords_4P_0_y, Coords_4P_0_z))
			myAssembly.features.changeKey(fromName='RP-1', toName='RP_Qz_4P_0')
			RP_Qz_4P_0 = (myAssembly.referencePoints.findAt((Coords_4P_0_x, Coords_4P_0_y, Coords_4P_0_z)),)
			RP_Qz_4P_0 = myAssembly.Set(referencePoints=RP_Qz_4P_0, name='RP_Qz_4P_0')
			
			# Load closer to x = L
			# Node numbers at x location of F
			nN_F_upchord 	= nNH*nNB*int(round((Pos_Qz_4P_L)/sizeEL,0))+int(nEB/2)+1				    # Node number at x location of F in the mid of the upper chord
			nN_F_mid 		= nNH*nNB*int(round((Pos_Qz_4P_L)/sizeEL,0))+int(nEH/2)*nNB+int(nEB/2)+1	# Node number at x location of F in the centroid
			nN_F_lowchord	= nNH*nNB*int(round((Pos_Qz_4P_L)/sizeEL,0))+nEH*nNB+int(nEB/2)+1		    # Node number at x location of F in the mid of the lower chord	
			rot_F = (myNodesCoords[nN_F_upchord-1][1] - myNodesCoords[nN_F_lowchord-1][1])/(myNodesCoords[nN_F_lowchord-1][2] - myNodesCoords[nN_F_upchord-1][2])	# Cross-section rotation at F due to imperfections
			
			# Coordinates of the reference point including effects of imperfections
			Coords_4P_L_x = myNodesCoords[nN_F_mid-1][0]
			Coords_4P_L_y = myNodesCoords[nN_F_mid-1][1] + ey_Qz * math.cos(rot_F) - ez_Qz * math.sin(rot_F)
			Coords_4P_L_z = myNodesCoords[nN_F_mid-1][2] + ez_Qz * math.cos(rot_F) + ey_Qz * math.sin(rot_F) 
			
			# Reference point
			myAssembly.ReferencePoint(point=(Coords_4P_L_x, Coords_4P_L_y, Coords_4P_L_z))
			myAssembly.features.changeKey(fromName='RP-1', toName='RP_Qz_4P_L')
			RP_Qz_4P_L = (myAssembly.referencePoints.findAt((Coords_4P_L_x, Coords_4P_L_y, Coords_4P_L_z)),)
			RP_Qz_4P_L = myAssembly.Set(referencePoints=RP_Qz_4P_L, name='RP_Qz_4P_L')
			
			if Moment == 'Qz_4P_area_T':
				Coords_4P_T_x = (Coords_4P_0_x + Coords_4P_L_x) / 2
				Coords_4P_T_y = (Coords_4P_0_y + Coords_4P_L_y) / 2
				Coords_4P_T_z = (Coords_4P_0_z + Coords_4P_L_z) / 2
				myAssembly.ReferencePoint(point=(Coords_4P_T_x, Coords_4P_T_y, Coords_4P_T_z))
				myAssembly.features.changeKey(fromName='RP-1', toName='RP_Qz_4P_area_T')
				RP_Qz_4P_area_T = (myAssembly.referencePoints.findAt((Coords_4P_T_x, Coords_4P_T_y, Coords_4P_T_z)),)
				RP_Qz_4P_area_T = myAssembly.Set(referencePoints=RP_Qz_4P_area_T, name='RP_Qz_4P_area_T')
			
		# 2 point loads - 4-point bending - distributed over an area
		if Moment == 'Qz_4P_area':
			# Load closer to x = 0
			# Node numbers at x location of F
			nN_F_upchord 	= nNH*nNB*int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0))+int(nEB/2)+1				    # Node number at x location of F in the mid of the upper chord
			nN_F_mid 		= nNH*nNB*int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0))+int(nEH/2)*nNB+int(nEB/2)+1	# Node number at x location of F in the centroid
			nN_F_lowchord	= nNH*nNB*int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0))+nEH*nNB+int(nEB/2)+1		    # Node number at x location of F in the mid of the lower chord	
			rot_F = (myNodesCoords[nN_F_upchord-1][1] - myNodesCoords[nN_F_lowchord-1][1])/(myNodesCoords[nN_F_lowchord-1][2] - myNodesCoords[nN_F_upchord-1][2])	# Cross-section rotation at F due to imperfections
			
			# Coordinates of the reference points including effects of imperfections
			Coords_4P_0_x = myNodesCoords[nN_F_mid-1][0]
			Coords_4P_0_y = myNodesCoords[nN_F_mid-1][1] + ey_Qz * math.cos(rot_F) - ez_Qz * math.sin(rot_F)
			Coords_4P_0_z = myNodesCoords[nN_F_mid-1][2] + ez_Qz * math.cos(rot_F) + ey_Qz * math.sin(rot_F)
					
			myAssembly.ReferencePoint(point=(Coords_4P_0_x, Coords_4P_0_y, Coords_4P_0_z))
			myAssembly.features.changeKey(fromName='RP-1', toName='RP_Qz_4P_0')
			RP_Qz_4P_0 = (myAssembly.referencePoints.findAt((Coords_4P_0_x, Coords_4P_0_y, Coords_4P_0_z)),)
			RP_Qz_4P_0 = myAssembly.Set(referencePoints=RP_Qz_4P_0, name='RP_Qz_4P_0')
			
			# Load closer to x = L
			# Node numbers at x location of F
			nN_F_upchord 	= nNH*nNB*int(round(Pos_Qz_4P_L/sizeEL,0))+int(nEB/2)+1				        # Node number at x location of F in the mid of the upper chord
			nN_F_mid 		= nNH*nNB*int(round(Pos_Qz_4P_L/sizeEL,0))+int(nEH/2)*nNB+int(nEB/2)+1		# Node number at x location of F in the centroid
			nN_F_lowchord	= nNH*nNB*int(round(Pos_Qz_4P_L/sizeEL,0))+nEH*nNB+int(nEB/2)+1		        # Node number at x location of F in the mid of the lower chord	
			rot_F = (myNodesCoords[nN_F_upchord-1][1] - myNodesCoords[nN_F_lowchord-1][1])/(myNodesCoords[nN_F_lowchord-1][2] - myNodesCoords[nN_F_upchord-1][2])	# Cross-section rotation at F due to imperfections
			
			# Coordinates of the reference points including effects of imperfections
			Coords_4P_L_x = myNodesCoords[nN_F_mid-1][0]
			Coords_4P_L_y = myNodesCoords[nN_F_mid-1][1] + ey_Qz * math.cos(rot_F) - ez_Qz * math.sin(rot_F)
			Coords_4P_L_z = myNodesCoords[nN_F_mid-1][2] + ez_Qz * math.cos(rot_F) + ey_Qz * math.sin(rot_F)
			
			myAssembly.ReferencePoint(point=(Coords_4P_L_x, Coords_4P_L_y, Coords_4P_L_z))
			myAssembly.features.changeKey(fromName='RP-1', toName='RP_Qz_4P_L')
			RP_Qz_4P_L = (myAssembly.referencePoints.findAt((Coords_4P_L_x, Coords_4P_L_y, Coords_4P_L_z)),)
			RP_Qz_4P_L = myAssembly.Set(referencePoints=RP_Qz_4P_L, name='RP_Qz_4P_L')
		
		# End rotations causing constant My
		if Moment == 'My_const' and Load_type_My == 'w' and Supp_type == 'real':
			# Loading at x = 0
			myAssembly.ReferencePoint(point=(0.0, 0.0, 0.0))
			myAssembly.features.changeKey(fromName='RP-1', toName='RP_My_const_0')
			RP_My_const_0 = (myAssembly.referencePoints.findAt((0.0, 0.0, 0.0)),)
			RP_My_const_0 = myAssembly.Set(referencePoints=RP_My_const_0, name='RP_My_const_0')
			# Loading at x = L
			myAssembly.ReferencePoint(point=(L, 0.0, (t_up2+t_u2)/2))
			myAssembly.features.changeKey(fromName='RP-1', toName='RP_My_const_L')
			RP_My_const_L = (myAssembly.referencePoints.findAt((L, 0.0, (t_up2+t_u2)/2)),)
			RP_My_const_L = myAssembly.Set(referencePoints=RP_My_const_L, name='RP_My_const_L')
			
	########################################
	''' Sets  '''
	########################################
	''' Supports '''
	# Support at x = L/2 in x direction - if displacement controlled actions in x direction are applied, this is not used
	# Supp_x - Node a cross-section centroid
	myNodes = []
	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+0.5*nEL*10000000+0.5*nEB*10000+nEH*0.5*10),)))
	N_Supp_x = myAssembly.Set(name='N_Supp_x', nodes=myNodes)

	# Supp_y for preventing LTB
	if type_var == 'auto' and noLTB == 1:
		myNodes = []
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+0.5*nEB*10000+intH*10),))))	for intL in range(1,nEL) for intH in range(1,nEH)	]
		N_Supp_noLTB = myAssembly.Set(name='N_Supp_noLTB', nodes=myNodes)
	
	# Supp_y for preventing in-plane buckling
	if type_var == 'auto' and noIPB == 1:
		myNodes = []
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+0.5*nEB*10000+intH*10),))))	for intL in range(1,nEL) for intH in range(1,nEH)	]
		N_Supp_noIPB = myAssembly.Set(name='N_Supp_noIPB', nodes=myNodes)
		
	# Ideal supports at cross-section centroid
	if Supp_type == 'ideal':

		# Displacement controlled loading and Nx or My_const and torsional actions
		if Step_n_Mx != 0 and ((Step_n_Nx != 0 and Load_type_Nx == 'w') or (Step_n_My != 0 and Moment == 'My_const' and Load_type_My == 'w')):
		
			# Supp_0_xyz
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intH*10+intB*10000),))))	for intH in range(nEH+1) for intB in range (nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intHR*10+intB*10000),))))	for intHR in range(nEH) for intB in range (nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intH*10+intBR*10000),))))	for intH in range(nEH+1) for intBR in range(nEB) 	]
			N_Supp_0_xyz = myAssembly.Set(name='N_Supp_0_xyz', nodes=myNodes)	
			
			# Supp_L_xyz
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+intH*10+intB*10000),))))		for intH in range(nEH+1) for intB in range (nEB+1) 	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEL*10000000+intHR*10+intB*10000),))))	for intHR in range(nEH) for intB in range (nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEL*10000000+intH*10+intBR*10000),))))	for intH in range(nEH+1) for intBR in range(nEB) 	]
			N_Supp_L_xyz = myAssembly.Set(name='N_Supp_L_xyz', nodes=myNodes)	
		
		# Displacement controlled loading and Nx or My_const
		elif (Step_n_Nx != 0 and Load_type_Nx == 'w') or (Step_n_My != 0 and Moment == 'My_const' and Load_type_My == 'w'):
		
			# Support at x = 0
			# Supp_0_xy - vertical line of nodes creating fork support
			myNodes = []
			if (nEH % 2) == 0:
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEB*0.5*10000+intHU*10),))))	for intHU in range(int(nEH/2))] 
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEB*0.5*10000+intHO*10),))))	for intHO in range(int(nEH/2)+1,nEH+1)]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEB*0.5*10000+intHR*10),))))	for intHR in range(nEH) ]
			else:
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEB*0.5*10000+intH*10),))))	for intH in range(nEH+1)] 
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEB*0.5*10000+intHRU*10),))))	for intHRU in range(int(nEH/2)) ]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEB*0.5*10000+intHRO*10),))))	for intHRO in range(int(nEH/2)+1,nEH) ]
			N_Supp_0_xy = myAssembly.Set(name='N_Supp_0_xy', nodes=myNodes)	
				
			# Supp_0_x
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010++intH*10+intBU*10000),))))	for intH in range(nEH+1) for intBU in range(int(nEB/2))				]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010++intH*10+intBO*10000),))))	for intH in range(nEH+1) for intBO in range ((int(nEB/2))+1, nEB+1) 	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015++intHR*10+intBU*10000),))))	for intHR in range(nEH) for intBU in range(int(nEB/2)) 				]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015++intHR*10+intBO*10000),))))	for intHR in range(nEH) for intBO in range ((int(nEB/2))+1, nEB+1) 	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010++intH*10+intBR*10000),))))	for intH in range(nEH+1) for intBR in range(nEB) ]
			N_Supp_0_x = myAssembly.Set(name='N_Supp_0_x', nodes=myNodes)	
			
			# Supp_0_xyz
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEB*0.5*10000+nEH*0.5*10),))))] 
			N_Supp_0_xyz = myAssembly.Set(name='N_Supp_0_xyz', nodes=myNodes)
			
			# Support at x = L
			# Supp_L_xy
			myNodes = []
			if (nEH % 2) == 0:
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+nEB*0.5*10000+intHU*10),))))	for intHU in range(int(nEH/2))] 
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+nEB*0.5*10000+intHO*10),))))	for intHO in range(int(nEH/2)+1,nEH+1)]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEL*10000000+nEB*0.5*10000+intHR*10),))))	for intHR in range(nEH) ]
			else:
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+nEB*0.5*10000+intH*10),))))	for intH in range(nEH+1)] 
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEL*10000000+nEB*0.5*10000+intHRU*10),))))	for intHRU in range(int(nEH/2)) ]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEL*10000000+nEB*0.5*10000+intHRO*10),))))	for intHRO in range(int(nEH/2)+1,nEH) ]
			N_Supp_L_xy = myAssembly.Set(name='N_Supp_L_xy', nodes=myNodes)

			# Supp_L_x
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+intH*10+intBU*10000),))))	for intH in range(nEH+1) for intBU in range(int(nEB/2))				]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+intH*10+intBO*10000),))))	for intH in range(nEH+1) for intBO in range ((int(nEB/2))+1, nEB+1) 	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEL*10000000+intHR*10+intBU*10000),))))	for intHR in range(nEH) for intBU in range(int(nEB/2)) 				]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEL*10000000+intHR*10+intBO*10000),))))	for intHR in range(nEH) for intBO in range ((int(nEB/2))+1, nEB+1) 	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEL*10000000+intH*10+intBR*10000),))))	for intH in range(nEH+1) for intBR in range(nEB) ]
			N_Supp_L_x = myAssembly.Set(name='N_Supp_L_x', nodes=myNodes)
			
			# Supp_0_xyz
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+nEB*0.5*10000+nEH*0.5*10),))))] 
			N_Supp_L_xyz = myAssembly.Set(name='N_Supp_L_xyz', nodes=myNodes)
			
		# Torsional actions
		elif Step_n_Mx != 0:
		
			# Supp_0_yz
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intH*10+intB*10000),))))	for intH in range(nEH+1) for intB in range(nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intHR*10+intB*10000),))))	for intHR in range(nEH) for intB in range(nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intH*10+intBR*10000),))))	for intH in range(nEH+1) for intBR in range(nEB)	]
			N_Supp_0_yz = myAssembly.Set(name='N_Supp_0_yz', nodes=myNodes)	
			
			# Supp_L_yz
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+intH*10+intB*10000),))))		for intH in range(nEH+1) for intB in range(nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEL*10000000+intHR*10+intB*10000),))))	for intHR in range(nEH) for intB in range(nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEL*10000000+intH*10+intBR*10000),))))	for intH in range(nEH+1) for intBR in range(nEB)	]
			N_Supp_L_yz = myAssembly.Set(name='N_Supp_L_yz', nodes=myNodes)	

		else:
			
			# Support at x = 0
			# Supp_0_yz - vertical line of nodes creating fork support
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+0.5*nEB*10000+intH*10),))))	for intH in range(nEH+1)] 
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+0.5*nEB*10000+intHR*10),))))	for intHR in range(nEH) ]
			N_Supp_0_yz = myAssembly.Set(name='N_Supp_0_yz', nodes=myNodes)	
				
			# Supp_0_z
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intH*10+intK*10000),))))		for intH in range(nEH+1) for intK in range(int(nEB/2)) ]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intH*10+intI*10000),))))		for intH in range(nEH+1) for intI in range ((int(nEB/2))+1, nEB+1) ]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intHR*10+intBU*10000),))))	for intHR in range(nEH) for intBU in range(int(nEB/2)) ]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intHR*10+intBO*10000),))))	for intHR in range(nEH)  for intBO in range ((int(nEB/2))+1, nEB+1) ]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intH*10+intBR*10000),))))		for intH in range(nEH+1) for intBR in range(nEB) ]
			N_Supp_0_z = myAssembly.Set(name='N_Supp_0_z', nodes=myNodes)	

			# Support at x = L
			# Supp_L_yz
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+0.5*nEB*10000+intH*10),)))) 	for intH in range(nEH+1)] 
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEL*10000000+0.5*nEB*10000+intHR*10),)))) for intHR in range(nEH) ]
			N_Supp_L_yz = myAssembly.Set(name='N_Supp_L_yz', nodes=myNodes)

			# Supp_L_z
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+intH*10+intK*10000),))))		for intH in range(nEH+1) for intK in range(int(nEB/2)) ]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+intH*10+intI*10000),))))		for intH in range(nEH+1) for intI in range ((int(nEB/2))+1, nEB+1) ]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEL*10000000+intHR*10+intBU*10000),))))	for intHR in range(nEH) for intBO in range(int(nEB/2)) ]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEL*10000000+intHR*10+intBO*10000),))))	for intHR in range(nEH)  for intBU in range ((int(nEB/2))+1, nEB+1) ]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEL*10000000+intH*10+intBR*10000),))))	for intH in range(nEH+1) for intBR in range(nEB) ]
			N_Supp_L_z = myAssembly.Set(name='N_Supp_L_z', nodes=myNodes)
				
	# Realistic supporting conditions at lower edge of the beam	
	elif Supp_type == 'real':

		# Support at x = 0
		# Position of fork support
		NL_y_0 = int(round((Pos_supp_0-sizeEL/4-0.00001)/sizeEL,0))			# in case of .75 nNL is rounded to the upper integer (by -0.00001), otherweise one gets unsymmetric support conditions for x = 0 & L 
		# Min & Max Node of the support area over the length nNL
		if (nEL_supp_0 % 2) == 0:		# check if even
			minNL_z_0 = max(int(round(Pos_supp_0/sizeEL-nEL_supp_0/2,0)),0)
			maxNL_z_0 = max(int(round(Pos_supp_0/sizeEL+nEL_supp_0/2,0)),minNL_z_0+1)
		else:
			if Pos_supp_0 > sizeEL*int(round(Pos_supp_0/sizeEL,0)):
				minNL_z_0 = max(int(round(Pos_supp_0/sizeEL,0)) - (nEL_supp_0-1)/2,0)
				maxNL_z_0 = max(int(round(Pos_supp_0/sizeEL,0)) + 1 + (nEL_supp_0-1)/2,minNL_z_0+1)
			else:
				minNL_z_0 = max(int(round(Pos_supp_0/sizeEL,0)) - 1 - (nEL_supp_0-1)/2,0)
				maxNL_z_0 = max(int(round(Pos_supp_0/sizeEL,0)) + (nEL_supp_0-1)/2,minNL_z_0+1)

		# Supp_0_y
		# Torsional actions
		if Step_n_Mx != 0:
			myNodes = []		
			if (round((Pos_supp_0-0.00001)/sizeEL*2,0) % 2) == 0:			# Check if mid nodes should be used
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intH*10+intB*10000+NL_y_0*10000000),)))) 	for intH in range(nEH) for intB in range (nEB+1) 	]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intH*10+intB*10000+NL_y_0*10000000),))))	for intH in range(nEH) for intB in range (nEB+1)	]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intH*10+intBR*10000+NL_y_0*10000000),))))	for intH in range(nEH) for intBR in range(nEB) 		]
			else:
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intH*10+intB*10000+NL_y_0*10000000),)))) 	for intH in range(nEH) for intB in range (nEB) 	]
			N_Supp_0_y = myAssembly.Set(name='N_Supp_0_y', nodes=myNodes)
		
		# Other actions
		else:
			myNodes = []
			if (round((Pos_supp_0-0.00001)/sizeEL*2,0) % 2) == 0:			# Check if mid nodes should be used
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intH*10+int(nEB/2)*10000+NL_y_0*10000000),))),
				myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intH*10+int(nEB/2)*10000+NL_y_0*10000000),))))
				for intH in range(nEH)]
			else:
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intH*10+int(nEB/2)*10000+NL_y_0*10000000),))))	for intH in range(nEH)]
			N_Supp_0_y = myAssembly.Set(name='N_Supp_0_y', nodes=myNodes)
		
		# Supp_0_z
		if ((Step_n_Nx != 0 and Load_type_Nx == 'w') or (Step_n_My != 0 and Moment == 'My_const' and Load_type_My == 'w')) and minNL_z_0 == 0:
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEH*10+intB*10000+intL*10000000),))))		for intB in range(nEB+1) for intL in range(minNL_z_0+1, maxNL_z_0+1) 	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEH*10+intBR*10000+intL*10000000),))))	for intBR in range(nEB)  for intL in range(minNL_z_0+1, maxNL_z_0+1) 	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+nEH*10+intB*10000+intLR*10000000),))))	for intB in range(nEB+1) for intLR in range(minNL_z_0, maxNL_z_0)		]
			N_Supp_0_z = myAssembly.Set(name='N_Supp_0_z', nodes=myNodes)	
		else:
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEH*10+intB*10000+intL*10000000),))))		for intB in range(nEB+1) for intL in range(minNL_z_0, NL_y_0)		]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEH*10+intB*10000+intL*10000000),))))		for intB in range(nEB+1) for intL in range(NL_y_0+1, maxNL_z_0+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEH*10+intBR*10000+intL*10000000),))))	for intBR in range(nEB)  for intL in range(minNL_z_0, NL_y_0)		]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEH*10+intBR*10000+intL*10000000),))))	for intBR in range(nEB)  for intL in range(NL_y_0+1, maxNL_z_0+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+nEH*10+intB*10000+intLR*10000000),))))	for intB in range(nEB+1) for intLR in range(minNL_z_0, maxNL_z_0)	]
			N_Supp_0_z = myAssembly.Set(name='N_Supp_0_z', nodes=myNodes)	
			# Supp_0_xz
			myNodes = []
			if (round((Pos_supp_0-0.00001)/sizeEL*2,0) % 2) == 0:			# Check if mid node should be omitted
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEH*10+intB*10000+NL_y_0*10000000),))))	for intB in range(int(nEB/2))			]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEH*10+intB*10000+NL_y_0*10000000),))))	for intB in range(int(nEB/2)+1,nEB+1) 	]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEH*10+intBR*10000+NL_y_0*10000000),))))	for intBR in range(nEB)  				]
			else:
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEH*10+intB*10000+NL_y_0*10000000),))))	for intB in range(nEB+1) 				]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEH*10+intBR*10000+NL_y_0*10000000),))))	for intBR in range(int(nEB/2))  		]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEH*10+intBR*10000+NL_y_0*10000000),))))	for intBR in range(int(nEB/2)+1,nEB)  	]
			N_Supp_0_xz = myAssembly.Set(name='N_Supp_0_xz', nodes=myNodes)
			# Supp_0_xyz
			myNodes = []
			if (round((Pos_supp_0-0.00001)/sizeEL*2,0) % 2) == 0:			# Check if mid node should be used
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEH*10+int(nEB/2)*10000+NL_y_0*10000000),)))) for intH in range(nEH)]
			else:
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+nEH*10+int(nEB/2)*10000+NL_y_0*10000000),))))	]
			N_Supp_0_xyz = myAssembly.Set(name='N_Supp_0_xyz', nodes=myNodes)
		
		del NL_y_0
		del minNL_z_0
		del maxNL_z_0
		
		# Support at x = L
		# Position of fork support
		NL_y_L = int(round((Pos_supp_L-sizeEL/4)/sizeEL,0)) 
		# Min & Max Node of the support area over the length nNL	
		if (nEL_supp_L % 2) == 0:		# check if even
			minNL_z_L = min(int(round(Pos_supp_L/sizeEL-nEL_supp_L/2,0)),nNL-2)
			maxNL_z_L = min(int(round(Pos_supp_L/sizeEL+nEL_supp_L/2,0)),nNL-1)
		else:
			if Pos_supp_L > sizeEL*int(round(Pos_supp_L/sizeEL,0)):
				minNL_z_L = min(int(round(Pos_supp_L/sizeEL,0)) - (nEL_supp_L-1)/2,nNL-2)
				maxNL_z_L = min(int(round(Pos_supp_L/sizeEL,0)) + 1 + (nEL_supp_L-1)/2,nNL-1)
			else:
				minNL_z_L = min(int(round(Pos_supp_L/sizeEL,0)) - 1 - (nEL_supp_L-1)/2,nNL-2)
				maxNL_z_L = min(int(round(Pos_supp_L/sizeEL,0)) + (nEL_supp_L-1)/2,nNL-1)

		# Supp_L_y
		# Torsional actions
		if Step_n_Mx != 0:
			myNodes = []
			if (round(Pos_supp_L/sizeEL*2,0) % 2) == 0:		# Check if mid nodes should be used
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intH*10+intB*10000+NL_y_L*10000000),))))	for intH in range(nEH)  for intB in range (nEB+1) 	]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intHR*10+intB*10000+NL_y_L*10000000),))))	for intHR in range(nEH) for intB in range (nEB+1) 	]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intH*10+intBR*10000+NL_y_L*10000000),))))	for intH in range(nEH)  for intBR in range(nEB) 	]
			else:
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intH*10+intB*10000+NL_y_L*10000000),))))	for intH in range(nEH)  for intB in range (nEB+1) 	]
			N_Supp_L_y = myAssembly.Set(name='N_Supp_L_y', nodes=myNodes)
		
		# Other actions
		else:
			myNodes = []
			if (round(Pos_supp_L/sizeEL*2,0) % 2) == 0:		# Check if mid nodes should be used
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intH*10+int(nEB/2)*10000+NL_y_L*10000000),))),
				myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intH*10+int(nEB/2)*10000+NL_y_L*10000000),))))
				for intH in range(nEH) for intHR in range(nEH)]
			else:
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intH*10+int(nEB/2)*10000+NL_y_L*10000000),))))	for intH in range(nEH) for intHR in range(nEH)]
			N_Supp_L_y = myAssembly.Set(name='N_Supp_L_y', nodes=myNodes)
		
		# Supp_L_z
		if ((Step_n_Nx != 0 and Load_type_Nx == 'w') or (Step_n_My != 0 and Moment == 'My_const' and Load_type_My == 'w')) and maxNL_z_L == nEL:
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEH*10+intB*10000+intL*10000000),))))		for intB in range(nEB+1) for intL in range(minNL_z_L, maxNL_z_L)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEH*10+intBR*10000+intL*10000000),))))	for intBR in range(nEB)  for intL in range(minNL_z_L, maxNL_z_L)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+nEH*10+intB*10000+intLR*10000000),))))	for intB in range(nEB+1) for intLR in range(minNL_z_L, maxNL_z_L)	]
			N_Supp_L_z = myAssembly.Set(name='N_Supp_L_z', nodes=myNodes)	
		else:
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEH*10+intB*10000+intL*10000000),))))		for intB in range(nEB+1) for intL in range(minNL_z_L, NL_y_L)]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEH*10+intB*10000+intL*10000000),))))		for intB in range(nEB+1) for intL in range(NL_y_L+1, maxNL_z_L+1)]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEH*10+intBR*10000+intL*10000000),))))	for intBR in range(nEB)  for intL in range(minNL_z_L, NL_y_L)]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEH*10+intBR*10000+intL*10000000),))))	for intBR in range(nEB)  for intL in range(NL_y_L+1, maxNL_z_L+1)]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+nEH*10+intB*10000+intLR*10000000),))))	for intB in range(nEB+1) for intLR in range(minNL_z_L, maxNL_z_L)	]
			N_Supp_L_z = myAssembly.Set(name='N_Supp_L_z', nodes=myNodes)		
			# Supp_L_xz
			myNodes = []
			if (round(Pos_supp_L/sizeEL*2,0) % 2) == 0:			# Check if mid node should be omitted
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEH*10+intB*10000+NL_y_L*10000000),))))	for intB in range(int(nEB/2))	]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEH*10+intB*10000+NL_y_L*10000000),))))	for intB in range(int(nEB/2)+1,nEB+1) 	]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEH*10+intBR*10000+NL_y_L*10000000),))))	for intBR in range(nEB)  	]
			else:
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEH*10+intB*10000+NL_y_L*10000000),))))	for intB in range(nEB+1) 	]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEH*10+intBR*10000+NL_y_L*10000000),))))	for intBR in range(int(nEB/2))  	]
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEH*10+intBR*10000+NL_y_L*10000000),))))	for intBR in range(int(nEB/2)+1,nEB)  	]
			N_Supp_L_xz = myAssembly.Set(name='N_Supp_L_xz', nodes=myNodes)
			# Supp_L_xyz
			myNodes = []
			if (round(Pos_supp_L/sizeEL*2,0) % 2) == 0:			# Check if mid node should be used
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEH*10+int(nEB/2)*10000+NL_y_L*10000000),)))) for intH in range(nEH)]
			else:
				[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+nEH*10+int(nEB/2)*10000+NL_y_L*10000000),))))	]
			N_Supp_L_xyz = myAssembly.Set(name='N_Supp_L_xyz', nodes=myNodes)
		
		del NL_y_L
		del minNL_z_L
		del maxNL_z_L
			
	''' Bracing '''
	# Supports of equivalent beam bracing
	if Brace_type == 'springs+beam':
		myNodes = []
		myNodes.append(myAssembly.instances['instEquBeamBracing'].nodes.sequenceFromLabels(labels=(int(1),)))
		N_Supp_brace_0 = myAssembly.Set(name='N_Supp_Brace_0', nodes=myNodes)
		myNodes = []
		myNodes.append(myAssembly.instances['instEquBeamBracing'].nodes.sequenceFromLabels(labels=(int(nEL+1),)))
		N_Supp_Brace_L = myAssembly.Set(name='N_Supp_Brace_L', nodes=myNodes)
		
	''' Actions '''
	# Axial actions
	if Step_n_Nx != 0:
		if Load_type_Nx == 'F':
			myElements = []
			[ myElements.append(myInstanceBeam.elements.sequenceFromLabels(labels=(int(10010010+intB*10000+intH*10), ))) for intB in range(nEB) for intH in range(nEH) ]
			myAssembly.Surface(face1Elements=myElements, name='SF_Nx_0')
			SF_Nx_0 = myAssembly.surfaces['SF_Nx_0']
			myElements = []
			[ myElements.append(myInstanceBeam.elements.sequenceFromLabels(labels=(int(10010010+(nEL-1)*10000000+intB*10000+intH*10), ))) for intB in range(nEB) for intH in range(nEH) ]
			myAssembly.Surface(face2Elements=myElements, name='SF_Nx_L')
			SF_Nx_L = myAssembly.surfaces['SF_Nx_L']	
	
		elif Load_type_Nx == 'Exp_LTB' and Supp_type == 'real':
			# Nx_0
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intH*10+intB*10000),))))	for intH in range(nEH/2-1,nEH/2+2) for intB in range (nEB+1) 	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intHR*10+intB*10000),))))for intHR in range(nEH/2-1,nEH/2+1)    for intB in range (nEB+1) 	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intH*10+intBR*10000),))))	for intH in range(nEH/2-1,nEH/2+2) for intBR in range(nEB) 	]
			N_Nx_0 = myAssembly.Set(name='N_Nx_0', nodes=myNodes)	
			
			# Nx_L
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+intH*10+intB*10000),))))	for intH in range(nEH/2-1,nEH/2+2) for intB in range (nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEL*10000000+intHR*10+intB*10000),))))	for intHR in range(nEH/2-1,nEH/2+1) for intB in range (nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEL*10000000+intH*10+intBR*10000),))))	for intH in range(nEH/2-1,nEH/2+2) for intBR in range(nEB) 	]
			N_Nx_L = myAssembly.Set(name='N_Nx_L', nodes=myNodes)
	
		elif Load_type_Nx == 'w' and Supp_type == 'real':
			# Nx_0
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intH*10+intB*10000),))))	for intH in range(nEH+1) for intB in range (nEB+1) 	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intHR*10+intB*10000),))))for intHR in range(nEH)    for intB in range (nEB+1) 	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intH*10+intBR*10000),))))	for intH in range(nEH+1) for intBR in range(nEB) 	]
			N_Nx_0 = myAssembly.Set(name='N_Nx_0', nodes=myNodes)	
			
			# Nx_L
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+intH*10+intB*10000),))))	for intH in range(nEH+1) for intB in range (nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEL*10000000+intHR*10+intB*10000),))))	for intHR in range(nEH) for intB in range (nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEL*10000000+intH*10+intBR*10000),))))	for intH in range(nEH+1) for intBR in range(nEB) 	]
			N_Nx_L = myAssembly.Set(name='N_Nx_L', nodes=myNodes)
	
	# Torsinal actions
	if Step_n_Mx != 0:
		myNodes = []
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Mx/sizeEL,0))*10000000+intH*10+intBl*10000),))))	for intH in range(nEH+1) 				for intBl in range(int(nEB/2)) 				]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Mx/sizeEL,0))*10000000+intH*10+intBo*10000),))))	for intH in range(nEH+1) 				for intBo in range ((int(nEB/2))+1, nEB+1) 	]		
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Mx/sizeEL,0))*10000000+intHl*10+int(nEB/2)*10000),))))	for intHl in range(int(nEH/2)) 														]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Mx/sizeEL,0))*10000000+intHo*10+int(nEB/2)*10000),))))	for intHo in range ((int(nEH/2))+1, nEH+1)											]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+int(round(Pos_Mx/sizeEL,0))*10000000+intHR*10+intB*10000),))))	for intHR in range(nEH) 				for intB in range(nEB+1) 				]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+int(round(Pos_Mx/sizeEL,0))*10000000+intH*10+intBR*10000),))))	for intH in range(nEH+1) 				for intBR in range(nEB) 				]
		N_Mx = myAssembly.Set(name='N_Mx', nodes=myNodes)
	
	# Actions causing My
	# Qz_3P
	if Step_n_My != 0 and Moment == 'Qz_3P':
		myNodes = []
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Mx/sizeEL,0))*10000000+intH*10+intBl*10000),))))	for intH in range(nEH+1) 				for intBl in range(int(nEB/2)) 				]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Mx/sizeEL,0))*10000000+intH*10+intBo*10000),))))	for intH in range(nEH+1) 				for intBo in range ((int(nEB/2))+1, nEB+1) 	]		
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Mx/sizeEL,0))*10000000+intHl*10+int(nEB/2)*10000),))))	for intHl in range(int(nEH/2)) 														]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Mx/sizeEL,0))*10000000+intHo*10+int(nEB/2)*10000),))))	for intHo in range ((int(nEH/2))+1, nEH+1)											]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+int(round(Pos_Mx/sizeEL,0))*10000000+intHR*10+intB*10000),))))	for intHR in range(nEH) 				for intB in range(nEB+1) 				]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+int(round(Pos_Mx/sizeEL,0))*10000000+intH*10+intBR*10000),))))	for intH in range(nEH+1) 				for intBR in range(nEB) 				]
		N_Qz_3P = myAssembly.Set(name='N_Qz_3P', nodes=myNodes)	
		
	# Qz_3Ps
	if Step_n_My != 0 and Moment == 'Qz_3Ps' and Load_type_My == 'F':
		myElements = []
		[ myElements.append(myInstanceBeam.elements.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Qz_3P/sizeEL-1,0))*10000000+intB*10000+intH*10), ))) for intB in range(nEB) for intH in range(nEH) ]
		myAssembly.Surface(face2Elements=myElements, name='SF_Qz_3Ps')
		SF_Qz_3Ps = myAssembly.surfaces['SF_Qz_3Ps']
		
		H_Qz_3Ps = lowchord[int(round(Pos_Qz_3P/sizeEL,0))]-upchord[int(round(Pos_Qz_3P/sizeEL,0))]	# Height of the surface SF_Qz_3Ps
		
	# Qz_3P_area
	if Step_n_My != 0 and Moment == 'Qz_3P_area':	
		if (nEL_area % 2) == 0:		# check if even
			min_nNL = int(round(Pos_Qz_3P/sizeEL,0) - nEL_area/2)
			max_nNL = int(round(Pos_Qz_3P/sizeEL,0) + nEL_area/2)
		else:
			if Pos_Qz_3P > sizeEL*int(round(Pos_Qz_3P/sizeEL,0)):
				min_nNL = int(round(Pos_Qz_3P/sizeEL,0) - (nEL_area-1)/2)
				max_nNL = int(round(Pos_Qz_3P/sizeEL,0) + 1 + (nEL_area-1)/2)
			else:
				min_nNL = int(round(Pos_Qz_3P/sizeEL,0) - 1 - (nEL_area-1)/2)
				max_nNL = int(round(Pos_Qz_3P/sizeEL,0) + (nEL_area-1)/2)
				
		myNodes = []
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Qz_3P/sizeEL,0))*10000000+int(0.5*nEB)*10000+0*10),))))]
		N_Qz_3P_area_yz = myAssembly.Set(name='N_Qz_3P_area_yz', nodes=myNodes)
		
		myNodes = []
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+int(0.5*nEB)*10000+0*10),))))	for intL in range(min_nNL,int(round(Pos_Qz_3P/sizeEL,0))) 	]
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+int(0.5*nEB)*10000+0*10),))))	for intL in range(int(round(Pos_Qz_3P/sizeEL,0))+1,max_nNL+1) 	]
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intLR*10000000+int(0.5*nEB)*10000+0*10),))))	for intLR in range(min_nNL,max_nNL)		]
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+0*10),))))	for intL in range(min_nNL,max_nNL+1) for intB in range(int(0.5*nEB)) 	]
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+0*10),))))	for intL in range(min_nNL,max_nNL+1) for intB in range(int(0.5*nEB)+1,nEB+1) 	]
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intLR*10000000+intB*10000+0*10),))))	for intLR in range(min_nNL,max_nNL)  for intB in range(int(0.5*nEB))	]
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intLR*10000000+intB*10000+0*10),))))	for intLR in range(min_nNL,max_nNL)  for intB in range(int(0.5*nEB)+1,nEB+1)	]
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intL*10000000+intJ*10000+0*10),))))	for intL in range(min_nNL,max_nNL+1) for intJ in range(nEB)		]
		N_Qz_3P_area_z = myAssembly.Set(name='N_Qz_3P_area_z', nodes=myNodes)
	
	# Exp_LTB
	if Step_n_My != 0 and Moment == 'Exp_LTB':	
		if (nEL_area % 2) == 0:		# check if even
			min_nNL = int(round(Pos_Qz_3P/sizeEL,0) - nEL_area/2)
			max_nNL = int(round(Pos_Qz_3P/sizeEL,0) + nEL_area/2)
		else:
			if Pos_Qz_3P > sizeEL*int(round(Pos_Qz_3P/sizeEL,0)):
				min_nNL = int(round(Pos_Qz_3P/sizeEL,0) - (nEL_area-1)/2)
				max_nNL = int(round(Pos_Qz_3P/sizeEL,0) + 1 + (nEL_area-1)/2)
			else:
				min_nNL = int(round(Pos_Qz_3P/sizeEL,0) - 1 - (nEL_area-1)/2)
				max_nNL = int(round(Pos_Qz_3P/sizeEL,0) + (nEL_area-1)/2)
				
		myNodes = []
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+0*10),))))	for intL in range(min_nNL,max_nNL+1) for intB in range(nEB+1) 	]
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intLR*10000000+intB*10000+0*10),))))	for intLR in range(min_nNL,max_nNL) 	 for intB in range(nEB+1)	]
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intL*10000000+intJ*10000+0*10),))))	for intL in range(min_nNL,max_nNL+1) for intJ in range(nEB)		]
		N_Qz_3P_area = myAssembly.Set(name='N_Qz_3P_area', nodes=myNodes)
		if Moment == 'Exp_LTB':
			myNodes = []
			[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intL*10000000+0*10000+intHR*10),))))	for intL in range(min_nNL,max_nNL+1) 	for intHR in range(1)	]
			N_Qz_3P_area_side = myAssembly.Set(name='N_Qz_3P_area_side', nodes=myNodes)
			myNodes = []
			for intL in range(min_nNL,max_nNL+1):
				for intB in range(1,nEB):
					for intH in range(1,max(2,int(L_screws_load/sizeEH)+1)):
						if int(10010010+intL*10000000+intB*10000+intH*10) != int(10010010+0.5*nEL*10000000+0.5*nEB*10000+nEH*0.5*10):
							myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+intH*10),)))
			for intLR in range(min_nNL,max_nNL):
				for intB in range(1,nEB):
					for intH in range(1,max(2,int(L_screws_load/sizeEH)+1)):
						if int(15010010+intLR*10000000+intB*10000+intH*10) != int(10010010+0.5*nEL*10000000+0.5*nEB*10000+nEH*0.5*10):
							myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intLR*10000000+intB*10000+intH*10),)))
			for intL in range(min_nNL,max_nNL+1):
				for intBR in range(nEB):
					for intH in range(1,max(2,int(L_screws_load/sizeEH)+1)):	
						if int(10015010+intL*10000000+intBR*10000+intH*10) != int(10010010+0.5*nEL*10000000+0.5*nEB*10000+nEH*0.5*10):
							myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intL*10000000+intBR*10000+intH*10),)))
			for intL in range(min_nNL,max_nNL+1):
				for intB in range(1,nEB):
					for intHR in range(0,max(1,int(L_screws_load/sizeEH))):
						if int(10010015+intL*10000000+intB*10000+intHR*10) != int(10010010+0.5*nEL*10000000+0.5*nEB*10000+nEH*0.5*10):
							myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intL*10000000+intB*10000+intHR*10),)))
			N_Qz_3P_screws = myAssembly.Set(name='N_Qz_3P_screws', nodes=myNodes)
			
	# Qz_4P
	if Step_n_My != 0 and Moment == 'Qz_4P' or Moment == 'Qz_4P_area_T':
		# Load closer to x = 0
		myNodes = []
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0))*10000000+intH*10+intBl*10000),))))	for intH in range(nEH+1)	for intBl in range(int(nEB/2))				]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0))*10000000+intH*10+intBo*10000),))))	for intH in range(nEH+1)	for intBo in range ((int(nEB/2))+1, nEB+1)	]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0))*10000000+intHl*10+int(nEB/2)*10000),))))	for intHl in range(int(nEH/2)) 											]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0))*10000000+intHo*10+int(nEB/2)*10000),))))	for intHo in range ((int(nEH/2))+1, nEH+1) 								]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0))*10000000+intHR*10+intB*10000),))))	for intHR in range(nEH)		for intB in range(nEB+1)				]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0))*10000000+intH*10+intBR*10000),))))	for intH in range(nEH+1)	for intBR in range(nEB)					]
		N_Qz_4P_0 = myAssembly.Set(name='N_Qz_4P_0', nodes=myNodes)	
		
		# Load closer to x = L
		myNodes = []
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Qz_4P_L/sizeEL,0))*10000000+intH*10+intBl*10000),))))	for intH in range(nEH+1)	for intBl in range(int(nEB/2))				]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Qz_4P_L/sizeEL,0))*10000000+intH*10+intBo*10000),))))	for intH in range(nEH+1)	for intBo in range ((int(nEB/2))+1, nEB+1)	]		
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Qz_4P_L/sizeEL,0))*10000000+intHl*10+int(nEB/2)*10000),))))	for intHl in range(int(nEH/2)) 											]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Qz_4P_L/sizeEL,0))*10000000+intHo*10+int(nEB/2)*10000),))))	for intHo in range ((int(int(nEH/2)))+1, nEH+1) 								]	
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+int(round(Pos_Qz_4P_L/sizeEL,0))*10000000+intHR*10+intB*10000),))))	for intHR in range(nEH)		for intB in range(nEB+1)				]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+int(round(Pos_Qz_4P_L/sizeEL,0))*10000000+intH*10+intBR*10000),))))	for intH in range(nEH+1)	for intBR in range(nEB)					]
		N_Qz_4P_L = myAssembly.Set(name='N_Qz_4P_L', nodes=myNodes)	
	
	# Qz_4Ps
	if Step_n_My != 0 and Moment == 'Qz_4Ps' and Load_type_My == 'F':
		# Load closer to x = 0
		myElements = []
		[ myElements.append(myInstanceBeam.elements.sequenceFromLabels(labels=(int(10010010+int(round((Pos_Qz_4P_0-0.00001)/sizeEL-1,0))*10000000+intB*10000+intH*10), ))) for intB in range(nEB) for intH in range(nEH) ]
		myAssembly.Surface(face2Elements=myElements, name='SF_Qz_4Ps_0')
		SF_Qz_4Ps_0 = myAssembly.surfaces['SF_Qz_4Ps_0']
		
		H_Qz_4Ps_0 = lowchord[int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0))]-upchord[int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0))]	# Height of the surface SF_Qz_4Ps_0
		
		# Load closer to x = L
		myElements = []
		[ myElements.append(myInstanceBeam.elements.sequenceFromLabels(labels=(int(10010010+int(round(Pos_Qz_4P_L/sizeEL-1,0))*10000000+intB*10000+intH*10), ))) for intB in range(nEB) for intH in range(nEH) ]
		myAssembly.Surface(face2Elements=myElements, name='SF_Qz_4Ps_L')
		SF_Qz_4Ps_L = myAssembly.surfaces['SF_Qz_4Ps_L']
		
		H_Qz_4Ps_L = lowchord[int(round(Pos_Qz_4P_L/sizeEL,0))]-upchord[int(round(Pos_Qz_4P_L/sizeEL,0))]	# Height of the surface SF_Qz_4Ps_L
		
	# Qz_4P_area
	if Step_n_My != 0 and (Moment == 'Qz_4P_area' or Moment == 'Qz_4P_area_T'):
		# Load closer to x = 0
		myNodes = []
		min_nNL_0 = int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0)) - max(1,round(nEL_area/2,0))
		max_nNL_0 = int(round((Pos_Qz_4P_0-0.00001)/sizeEL,0)) + max(1,round(nEL_area/2,0))
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+0*10),))))	for intL in range(min_nNL_0,max_nNL_0+1)	for intB in range(nEB+1)	]
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intLR*10000000+intB*10000+0*10),))))	for intLR in range(min_nNL_0,max_nNL_0) 	for intB in range(nEB+1) 	]
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intL*10000000+intJ*10000+0*10),))))	for intL in range(min_nNL_0,max_nNL_0+1) 	for intJ in range(nEB)		]
		N_Qz_4P_area_0 = myAssembly.Set(name='N_Qz_4P_area_0', nodes=myNodes)
		
		# Load closer to x = L
		myNodes = []
		min_nNL_L = int(round(Pos_Qz_4P_L/sizeEL,0)) - max(1,round(nEL_area/2,0))
		max_nNL_L = int(round(Pos_Qz_4P_L/sizeEL,0)) + max(1,round(nEL_area/2,0))
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+0*10),))))	for intL in range(min_nNL_L,max_nNL_L+1) 	for intB in range(nEB+1)	]
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intLR*10000000+intB*10000+0*10),))))	for intLR in range(min_nNL_L,max_nNL_L) 	for intB in range(nEB+1)	]
		[ (	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intL*10000000+intJ*10000+0*10),))))	for intL in range(min_nNL_L,max_nNL_L+1) 	for intJ in range(nEB)		]
		N_Qz_4P_area_L = myAssembly.Set(name='N_Qz_4P_area_L', nodes=myNodes)
	
	# My_const
	if Step_n_My != 0 and Moment == 'My_const':
		if Load_type_My == 'F':
			myElements = []
			[ myElements.append(myInstanceBeam.elements.sequenceFromLabels(labels=(int(10010010+intB*10000+intH*10), ))) for intB in range(nEB) for intH in range(nEH) ]
			myAssembly.Surface(face1Elements=myElements, name='SF_My_const_0')
			SF_My_const_0 = myAssembly.surfaces['SF_My_const_0']
			myElements = []
			[ myElements.append(myInstanceBeam.elements.sequenceFromLabels(labels=(int(10010010+(nEL-1)*10000000+intB*10000+intH*10), ))) for intB in range(nEB) for intH in range(nEH) ]
			myAssembly.Surface(face2Elements=myElements, name='SF_My_const_L')
			SF_My_const_L = myAssembly.surfaces['SF_My_const_L']	
			
		if Load_type_My == 'w' and Supp_type == 'real':
			# My_const_0
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intH*10+intB*10000),))))	for intH in range(nEH+1)	for intB in range (nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intHR*10+intB*10000),))))	for intHR in range(nEH)		for intB in range (nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intH*10+intM*10000),))))	for intH in range(nEH+1) 	for intM in range(nEB) 		]
			N_My_const_0 = myAssembly.Set(name='N_My_const_0', nodes=myNodes)	
		
			# My_const_L
			myNodes = []
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+intH*10+intB*10000),))))		for intH in range(nEH+1) 	for intB in range (nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEL*10000000+intHR*10+intB*10000),))))	for intHR in range(nEH) 	for intB in range (nEB+1)	]
			[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEL*10000000+intH*10+intM*10000),))))		for intH in range(nEH+1) 	for intM in range(nEB) 		]
			N_My_const_L = myAssembly.Set(name='N_My_const_L', nodes=myNodes)
		
	# qz
	if Step_n_My != 0 and Moment == 'qz' and Load_type_My == 'F':
		myElements = []
		[ myElements.append(myInstanceBeam.elements.sequenceFromLabels(labels=(int(10010010+intB*10000+intL*10000000), ))) for intB in range(nEB) for intL in range(nEL) ]
		myAssembly.Surface(face3Elements=myElements, name='SF_qz')
		SF_qz = myAssembly.surfaces['SF_qz']
		
	''' Output '''
	# Midspan - centroid
	myNodes = []
	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+0.5*nEL*10000000+0.5*nEB*10000+0.5*nEH*10),)))
	Out_Midspan_cen = myAssembly.Set(name='Out_Midspan_cen', nodes=myNodes)
	
	# Midspan - upper edge - centroid
	myNodes = []
	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+0.5*nEL*10000000+0.5*nEB*10000+0*10),)))
	Out_Midspan_upchord_cen = myAssembly.Set(name='Out_Midspan_upchord_cen', nodes=myNodes)
	
	# Midspan - lower edge - centroid
	myNodes = []
	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+0.5*nEL*10000000+0.5*nEB*10000+nEH*10),)))
	Out_Midspan_lowchord_cen = myAssembly.Set(name='Out_Midspan_lowchord_cen', nodes=myNodes)
	
	# Midspan - 1st cross-section edge of compression zone
	myNodes = []
	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+0.5*nEL*10000000+0*10000+0*10),)))
	Out_Midspan_upchord1 = myAssembly.Set(name='Out_Midspan_upchord1', nodes=myNodes)
	
	# Midspan - 2nd cross-section edge of compression zone
	myNodes = []
	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+0.5*nEL*10000000+nEB*10000+0*10),)))
	Out_Midspan_upchord2 = myAssembly.Set(name='Out_Midspan_upchord2', nodes=myNodes)
	
	# Midspan - 1st cross-section edge of tension zone
	myNodes = []
	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+0.5*nEL*10000000+0*10000+nEH*10),)))
	Out_Midspan_lowchord1 = myAssembly.Set(name='Out_Midspan_lowchord1', nodes=myNodes)
	
	# Midspan - 2nd cross-section edge of tension zone
	myNodes = []
	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+0.5*nEL*10000000+nEB*10000+nEH*10),)))
	Out_Midspan_lowchord2 = myAssembly.Set(name='Out_Midspan_lowchord2', nodes=myNodes)
	
	# Midspan - Side 3 centroid
	myNodes = []
	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+0.5*nEL*10000000+0*10000+0.5*nEH*10),)))
	Out_Midspan_cen_S3 = myAssembly.Set(name='Out_Midspan_cen_S3', nodes=myNodes)
	
	# Midspan - whole cross-section -> for M_y/z
	myNodes = []
	[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+0.5*nEL*10000000+intH*10+intB*10000),))))		for intH in range(nEH+1) 	for intB in range(nEB+1)	]
	[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+0.5*nEL*10000000+intHR*10+intB*10000),))))	for intHR in range (nEH)	for intB in range(nEB+1)	]
	[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+0.5*nEL*10000000+intH*10+intBR*10000),)))) 	for intH in range(nEH+1) 	for intBR in range(nEB)		]
	Out_Midspan_cs = myAssembly.Set(name='Out_Midspan_cs', nodes=myNodes)	
	
	# Midspan - middle of cross-section over the height -> for e.g. warping
	myNodes = []
	[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+0.5*nEL*10000000+intH*10+0.5*nEB*10000),))))	for intH in range(nEH+1)	]
	[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+0.5*nEL*10000000+intHR*10+0.5*nEB*10000),)))) 	for intHR in range (nEH)	]
	Out_Midspan_LineH = myAssembly.Set(name='Out_Midspan_LineH', nodes=myNodes)	
	
	# Midspan - 2.5H - centroid (for determination of local bending stiffness)
	if Moment == 'Qz_4P_area_T':
		myNodes = []
		myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+(0.5*nEL-2.5*H/sizeEL)*10000000+0.5*nEB*10000+0.5*nEH*10),)))
		Out_Eloc_1_cen = myAssembly.Set(name='Out_Eloc_1_cen', nodes=myNodes)
		myNodes = []
		myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+(0.5*nEL+2.5*H/sizeEL)*10000000+0.5*nEB*10000+0.5*nEH*10),)))
		Out_Eloc_2_cen = myAssembly.Set(name='Out_Eloc_2_cen', nodes=myNodes)
	
	# Support_0 - whole cross-section -> for M_y/z
	myNodes = []
	[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+0*nEL*10000000+intH*10+intB*10000),))))		for intH in range(nEH+1) 	for intB in range(nEB+1)	]
	[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+0*nEL*10000000+intHR*10+intB*10000),))))	for intHR in range (nEH)	for intB in range(nEB+1)	]
	[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+0*nEL*10000000+intH*10+intBR*10000),))))	for intH in range(nEH+1) 	for intBR in range(nEB)		]
	Out_Supp_0_cs = myAssembly.Set(name='Out_Supp_0_cs', nodes=myNodes)
	
	# Support_L - whole cross-section -> for M_y/z
	myNodes = []
	[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+nEL*10000000+intH*10+intB*10000),))))	for intH in range(nEH+1) 	for intB in range(nEB+1)	]
	[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+nEL*10000000+intHR*10+intB*10000),))))	for intHR in range (nEH)	for intB in range(nEB+1)	]
	[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+nEL*10000000+intH*10+intBR*10000),))))	for intH in range(nEH+1) 	for intBR in range(nEB)		]
	Out_Supp_L_cs = myAssembly.Set(name='Out_Supp_L_cs', nodes=myNodes)
	
	# Tension Zone
	myElements = []
	[	myElements.append(myAssembly.instances['instbeam'].elements.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+intH*10),)))	for intL in range(nEL) for intB in range(nEB) for intH in range(int(nEH/2),nEH)	]
	Out_TensZone_all = myAssembly.Set(name='Out_TensZone_all', elements=myElements)
	
	# Nodes for tensile failure
	if elem_deletion == 0:
		myNodes = []
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+intH*10+0*10000),))))		for intL in range(int(nEL/4),int(nNL-nEL/4)) 	for intH in range(nEH+1) 	]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+intH*10+nEB*10000),))))		for intL in range(int(nEL/4),int(nNL-nEL/4)) 	for intH in range(nEH+1) 	]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intLR*10000000+intH*10+0*10000),))))		for intLR in range(int(nEL/4),int(nNL-nEL/4-1)) for intH in range(nEH+1) 	]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intLR*10000000+intH*10+nEB*10000),))))	for intLR in range(int(nEL/4),int(nNL-nEL/4-1)) for intH in range(nEH+1)	]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intL*10000000+intHR*10+0*10000),))))		for intL in range(int(nEL/4),int(nNL-nEL/4)) 	for intHR in range(nEH)		]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+intL*10000000+intHR*10+nEB*10000),))))	for intL in range(int(nEL/4),int(nNL-nEL/4)) 	for intHR in range(nEH)		]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+0*10+intB*10000),))))		for intL in range(int(nEL/4),int(nNL-nEL/4)) 	for intB in range(1,nEB) 	]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+nEH*10+intB*10000),))))		for intL in range(int(nEL/4),int(nNL-nEL/4)) 	for intB in range(1,nEB) 	]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intLR*10000000+0*10+intB*10000),))))		for intLR in range(int(nEL/4),int(nNL-nEL/4-1)) for intB in range(1,nEB) 	]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(15010010+intLR*10000000+nEH*10+intB*10000),))))	for intLR in range(int(nEL/4),int(nNL-nEL/4-1)) for intB in range(1,nEB)	]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intL*10000000+0*10+intBR*10000),))))		for intL in range(int(nEL/4),int(nNL-nEL/4)) 	for intBR in range(nEB)		]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+intL*10000000+nEH*10+intBR*10000),))))	for intL in range(int(nEL/4),int(nNL-nEL/4)) 	for intBR in range(nEB)		]
		Out_TensZone_fail_N = myAssembly.Set(name='Out_TensZone_fail_N', nodes=myNodes)
	# Elements for tensile failure
	else:
		myElements = []
		[	myElements.append(myAssembly.instances['instbeam'].elements.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+intH*10),)))	for intL in range(int(nEL/8),int(nEL-nEL/8)) for intB in range(nEB) for intH in range(nEH-1,nEH)	]
		Out_TensZone_fail_E = myAssembly.Set(name='Out_TensZone_fail_E', elements=myElements)

	# Nodes for compression failure
	myNodes = []
	[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+intH*10+intB*10000),))))	for intL in range(int(nEL/8),int(nNL-nEL/8)) for intH in range(nEH+1) 	for intB in range(nEB+1)	]
	Out_CompZone_fail_N = myAssembly.Set(name='Out_CompZone_fail_N', nodes=myNodes)
		
	# # Nodes for shear failure
	nEH_red = int(max(round(nEH/6.,0),1))
	if elem_deletion == 0:
		myNodes = []
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+intH*10),))))	for intL in range(int(nEL_supp_shear_0),int(nEL/2-nEL_screws_load/2+1)) 		for intH in range(nEH_red,nEH+1-nEH_red) 	for intB in range(1) 	]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+intH*10),))))	for intL in range(int(nEL/2+nEL_screws_load/2),int(nEL-nEL_supp_shear_L+1)) 	for intH in range(nEH_red,nEH+1-nEH_red) 	for intB in range(1)	]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+intH*10),))))	for intL in range(int(nEL_supp_shear_0),int(nEL/2-nEL_screws_load/2+1)) 		for intH in range(nEH_red,nEH+1-nEH_red) 	for intB in range(nEB,nEB+1) 	]
		[ ( myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+intH*10),))))	for intL in range(int(nEL/2+nEL_screws_load/2),int(nEL-nEL_supp_shear_L+1)) 	for intH in range(nEH_red,nEH+1-nEH_red) 	for intB in range(nEB,nEB+1)	]
		Out_ShearZone_fail_N = myAssembly.Set(name='Out_ShearZone_fail_N', nodes=myNodes)
	# Elements for shear failure
	else:
		myElements = []
		[ ( myElements.append(myAssembly.instances['instbeam'].elements.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+intH*10),))))	for intL in range(int(nEL_supp_shear_0),int(nEL/2-nEL_screws_load/2)) 		for intH in range(nEH_red,nEH-nEH_red) 	for intB in range(1) 	]
		[ ( myElements.append(myAssembly.instances['instbeam'].elements.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+intH*10),))))	for intL in range(int(nEL/2+nEL_screws_load/2),int(nEL-nEL_supp_shear_L)) 	for intH in range(nEH_red,nEH-nEH_red) 	for intB in range(1)	]
		[ ( myElements.append(myAssembly.instances['instbeam'].elements.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+intH*10),))))	for intL in range(int(nEL_supp_shear_0),int(nEL/2-nEL_screws_load/2)) 		for intH in range(nEH_red,nEH-nEH_red) 	for intB in range(nEB-1,nEB) 	]
		[ ( myElements.append(myAssembly.instances['instbeam'].elements.sequenceFromLabels(labels=(int(10010010+intL*10000000+intB*10000+intH*10),))))	for intL in range(int(nEL/2+nEL_screws_load/2),int(nEL-nEL_supp_shear_L)) 	for intH in range(nEH_red,nEH-nEH_red) 	for intB in range(nEB-1,nEB)	]
		Out_ShearZone_fail_E = myAssembly.Set(name='Out_ShearZone_fail_E', elements=myElements)
	del nEH_red
		
	########################################
	''' Boundary Conditions '''
	########################################
	''' Supports '''
	# x Support
	myModel.DisplacementBC(name='Supp_x', createStepName='Initial', region=N_Supp_x, u1=SET, u2=UNSET, u3=UNSET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

	# y Support for preventing LTB
	if type_var == 'auto' and noLTB == 1 and noIPB == 1:
		myModel.DisplacementBC(name='Supp_noLTB', createStepName='Initial', region=N_Supp_noLTB, u1=UNSET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
	elif type_var == 'auto' and noLTB == 1 and noIPB == 0:
		myModel.DisplacementBC(name='Supp_noLTB', createStepName='Initial', region=N_Supp_noLTB, u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
	elif type_var == 'auto' and noLTB == 0 and noIPB == 1:
		myModel.DisplacementBC(name='Supp_noIPB', createStepName='Initial', region=N_Supp_noIPB, u1=UNSET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

	# z Support
	if Supp_type == 'ideal':
		if (Step_n_Nx != 0 and Load_type_Nx == 'w') or (Step_n_My != 0 and Moment == 'My_const' and Load_type_My == 'w'):
			myModel.DisplacementBC(name='Supp_0', createStepName='Initial', region=RP_Supp_0, u1=UNSET, u2=SET, u3=SET, ur1=SET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			myModel.DisplacementBC(name='Supp_L', createStepName='Initial', region=RP_Supp_L, u1=UNSET, u2=SET, u3=SET, ur1=SET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
		else:
			myModel.DisplacementBC(name='Supp_0', createStepName='Initial', region=RP_Supp_0, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			myModel.DisplacementBC(name='Supp_L', createStepName='Initial', region=RP_Supp_L, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
	else:
		myModel.DisplacementBC(name='Supp_0', createStepName='Initial', region=RP_Supp_0, u1=UNSET, u2=SET, u3=SET, ur1=SET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
		myModel.DisplacementBC(name='Supp_L', createStepName='Initial', region=RP_Supp_L, u1=UNSET, u2=SET, u3=SET, ur1=SET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
		
	if Supp_type == 'real':
		# Support at x = 0
		myModel.DisplacementBC(name='ForkSupp_0', createStepName='Initial', region=RP_ForkSupp_0, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
		# Support at x = L
		myModel.DisplacementBC(name='ForkSupp_L', createStepName='Initial', region=RP_ForkSupp_L, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
	
	''' Bracing '''
	if Brace_type == 'springs+beam':
		myModel.DisplacementBC(name='Supp_BeamBrace_0', createStepName='Initial', region=myModel.rootAssembly.sets['N_Supp_Brace_0'], u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
		myModel.DisplacementBC(name='Supp_BeamBrace_L', createStepName='Initial', region=myModel.rootAssembly.sets['N_Supp_Brace_L'], u1=UNSET, u2=UNSET, u3=SET, ur1=SET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
		
		# DOF of the bracings' RPs which are not supported by the springs are fixed
		for intA in range(len(pos_brace)):		
			myNameSupp = 'Supp_Brace_'+str(intA)
			if intA == 0:
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_0, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 1:
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_1, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 2:
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_2, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 3:
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_3, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 4:
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_4, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 5:
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_5, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 6:
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_6, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 7:
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_7, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 8:
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_8, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 9:
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_9, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

	'''	Actions '''
	# Torsional actions
	if Step_n_Mx != 0:
		myModel.DisplacementBC(name='Supp_Mx', createStepName='Initial', region=RP_Mx, u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
	
	# Actions causing My
	if Step_n_My != 0:
		# Point load - 3-point bending
		if Moment == 'Qz_3P': 
			myModel.DisplacementBC(name='Supp_Qz_3P', createStepName='Initial', region=RP_Qz_3P, u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
		
		# Point load - 3-point bending - distributed over an area
		if Moment == 'Qz_3P_area' or Moment == 'Exp_LTB': 
			myModel.DisplacementBC(name='Supp_Qz_3P_area', createStepName='Initial', region=RP_Qz_3P_area, u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
		
		# Point load - 4-point bending		
		if Moment == 'Qz_4P': 
			myModel.DisplacementBC(name='Supp_Qz_4P_0', createStepName='Initial', region=RP_Qz_4P_0, u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			myModel.DisplacementBC(name='Supp_Qz_4P_L', createStepName='Initial', region=RP_Qz_4P_L, u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
		
		# Point load - 4-point bending - distributed over an area
		if Moment == 'Qz_4P_area': 
			myModel.DisplacementBC(name='Supp_Qz_4P_area_0', createStepName='Initial', region=RP_Qz_4P_0, u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			myModel.DisplacementBC(name='Supp_Qz_4P_area_L', createStepName='Initial', region=RP_Qz_4P_L, u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
		if Moment == 'Qz_4P_area_T': 
			myModel.DisplacementBC(name='Supp_Qz_4P_area_T', createStepName='Initial', region=RP_Qz_4P_area_T, u1=SET, u2=UNSET, u3=UNSET, ur1=SET, ur2=UNSET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
	
	########################################
	''' Fricton at supports '''
	########################################
	if Step_n_Nx != 0 and M_mu_supp != 0:
		# Connector properties
		mdb.models[modelName].ConnectorSection(name='ConSec_Spring', rotationalType=ROTATION)
		mdb.models[modelName].sections['ConSec_Spring'].setValues(behaviorOptions =(connectorBehavior.ConnectorElasticity(components=(5,6, ), behavior=NONLINEAR, table=((-M_mu_supp, -0.0001), (M_mu_supp, 0.0001))), ) )

		# Feature
		if Supp_type == 'ideal':
			wire = myAssembly.WirePolyLine(points=((myAssembly.referencePoints[3], None), ), mergeType=IMPRINT, meshable=False)
		else:
			wire = myAssembly.WirePolyLine(points=((myAssembly.referencePoints[11], None), ), mergeType=IMPRINT, meshable=False)
		oldName = wire.name
		myAssembly.features.changeKey(fromName=oldName, toName='Wire_Spring_0')
		if Supp_type == 'ideal':
			wire = myAssembly.WirePolyLine(points=((myAssembly.referencePoints[5], None), ), mergeType=IMPRINT, meshable=False)
		else:
			wire = myAssembly.WirePolyLine(points=((myAssembly.referencePoints[13], None), ), mergeType=IMPRINT, meshable=False)
		oldName = wire.name
		myAssembly.features.changeKey(fromName=oldName, toName='Wire_Spring_L')

		# Assign Feature 'Wire' to RP
		myAssembly.Set(edges=myAssembly.edges.findAt(((0.0, y_Nx_0, z_Supp_ideal_0+ez_Nx_0),),), name='Wire_Spring_0')
		myAssembly.Set(edges=myAssembly.edges.findAt(((L, ey_Nx_L, z_Supp_ideal_L+ez_Nx_L),),), name='Wire_Spring_L')

		# Assign connector properties to 'Wire'
		Wire_Spring_0 = myAssembly.SectionAssignment(sectionName='ConSec_Spring', region=myAssembly.sets['Wire_Spring_0'])
		Wire_Spring_L = myAssembly.SectionAssignment(sectionName='ConSec_Spring', region=myAssembly.sets['Wire_Spring_L'])
	
	########################################
	''' Springs of the bracing '''
	########################################
	if Brace_type == 'springs+beam':
		# Eaves purlin at x = 0
		myNodes = myAssembly.instances['instEquBeamBracing'].nodes.sequenceFromLabels(labels=(int(1),))
		myName = 'SP_EPurlin_0'
		myAssembly.Set(name=myName, nodes=myNodes)
		region=((myAssembly.sets[myName], myAssembly.sets['RP_Supp_0']), )
		myAssembly.engineeringFeatures.TwoPointSpringDashpot(name=myName, regionPairs=region, axis=FIXED_DOF, dof1=2, dof2=2, orientation=None, springBehavior=ON, springStiffness=K_epurlin_0, dashpotBehavior=OFF, dashpotCoefficient=0.0)

		# Eaves purlin at x = L
		myNodes = myAssembly.instances['instEquBeamBracing'].nodes.sequenceFromLabels(labels=(int(nNL),))
		myName = 'SP_EPurlin_L'
		myAssembly.Set(name=myName, nodes=myNodes)
		region=((myAssembly.sets[myName], myAssembly.sets['RP_Supp_L']), )
		myAssembly.engineeringFeatures.TwoPointSpringDashpot(name=myName, regionPairs=region, axis=FIXED_DOF, dof1=2, dof2=2, orientation=None, springBehavior=ON, springStiffness=K_epurlin_0, dashpotBehavior=OFF, dashpotCoefficient=0.0)
	
		# Bracing
		for intA in range(len(pos_brace)):	
			myNodes = myAssembly.instances['instEquBeamBracing'].nodes.sequenceFromLabels(labels=(int(round(nEL*pos_brace[intA],0)+1),))
			myNameBeam = 'SP_Brace_'+str(intA)
			myAssembly.Set(name=myNameBeam, nodes=myNodes)			
			myName = 'RP_Brace_'+str(intA)
			region=((myAssembly.sets[myName], myAssembly.sets[myNameBeam]), )
			myAssembly.engineeringFeatures.TwoPointSpringDashpot(name=myNameBeam, regionPairs=region, axis=FIXED_DOF, dof1=2, dof2=2, orientation=None, springBehavior=ON, springStiffness=K_beam_brace[intA], dashpotBehavior=OFF, dashpotCoefficient=0.0)
			
	elif Brace_type == 'NodestoGround':
		for intA in range(len(pos_brace)):					
			myName = 'Brace_'+str(intA)
			myNameSupp = 'Supp_Brace_'+str(intA)
			if intA == 0:
				myAssembly.engineeringFeatures.SpringDashpotToGround(name=myName, region=RP_Brace_0, orientation=None, dof=2, springBehavior=ON, springStiffness=K_NodestoGround[intA], dashpotBehavior=OFF, dashpotCoefficient=0.0)							
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_0, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 1:
				myAssembly.engineeringFeatures.SpringDashpotToGround(name=myName, region=RP_Brace_1, orientation=None, dof=2, springBehavior=ON, springStiffness=K_NodestoGround[intA], dashpotBehavior=OFF, dashpotCoefficient=0.0)
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_1, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 2:
				myAssembly.engineeringFeatures.SpringDashpotToGround(name=myName, region=RP_Brace_2, orientation=None, dof=2, springBehavior=ON, springStiffness=K_NodestoGround[intA], dashpotBehavior=OFF, dashpotCoefficient=0.0)
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_2, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 3:
				myAssembly.engineeringFeatures.SpringDashpotToGround(name=myName, region=RP_Brace_3, orientation=None, dof=2, springBehavior=ON, springStiffness=K_NodestoGround[intA], dashpotBehavior=OFF, dashpotCoefficient=0.0)
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_3, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 4:
				myAssembly.engineeringFeatures.SpringDashpotToGround(name=myName, region=RP_Brace_4, orientation=None, dof=2, springBehavior=ON, springStiffness=K_NodestoGround[intA], dashpotBehavior=OFF, dashpotCoefficient=0.0)
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_4, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 5:
				myAssembly.engineeringFeatures.SpringDashpotToGround(name=myName, region=RP_Brace_5, orientation=None, dof=2, springBehavior=ON, springStiffness=K_NodestoGround[intA], dashpotBehavior=OFF, dashpotCoefficient=0.0)
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_5, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 6:
				myAssembly.engineeringFeatures.SpringDashpotToGround(name=myName, region=RP_Brace_6, orientation=None, dof=2, springBehavior=ON, springStiffness=K_NodestoGround[intA], dashpotBehavior=OFF, dashpotCoefficient=0.0)
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_6, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 7:
				myAssembly.engineeringFeatures.SpringDashpotToGround(name=myName, region=RP_Brace_7, orientation=None, dof=2, springBehavior=ON, springStiffness=K_NodestoGround[intA], dashpotBehavior=OFF, dashpotCoefficient=0.0)
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_7, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 8:
				myAssembly.engineeringFeatures.SpringDashpotToGround(name=myName, region=RP_Brace_8, orientation=None, dof=2, springBehavior=ON, springStiffness=K_NodestoGround[intA], dashpotBehavior=OFF, dashpotCoefficient=0.0)
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_8, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			elif intA == 9:
				myAssembly.engineeringFeatures.SpringDashpotToGround(name=myName, region=RP_Brace_9, orientation=None, dof=2, springBehavior=ON, springStiffness=K_NodestoGround[intA], dashpotBehavior=OFF, dashpotCoefficient=0.0)
				myModel.DisplacementBC(name=myNameSupp, createStepName='Initial', region=RP_Brace_9, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)		
	
	########################################
	''' Couplings '''
	########################################
	''' Supports '''
	# Ideal supports at cross-section centroid
	if Supp_type == 'ideal':
		
		# Displacement controlled loading and Nx or My_const and torsional actions
		if Step_n_Mx != 0 and ((Step_n_Nx != 0 and Load_type_Nx == 'w') or (Step_n_My != 0 and Moment == 'My_const' and Load_type_My == 'w')):
			myModel.Coupling(name='CP_Supp_0_xyz', controlPoint=RP_Supp_0, surface=N_Supp_0_xyz, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=ON)
			myModel.Coupling(name='CP_Supp_L_xyz', controlPoint=RP_Supp_L, surface=N_Supp_L_xyz, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=ON)
		
		# Displacement controlled loading and Nx or My_const
		elif (Step_n_Nx != 0 and Load_type_Nx == 'w') or (Step_n_My != 0 and Moment == 'My_const' and Load_type_My == 'w'):
			myModel.Coupling(name='CP_Supp_0_x', controlPoint=RP_Supp_0, surface=N_Supp_0_x, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=OFF, u3=OFF, ur1=OFF, ur2=OFF, ur3=ON)
			myModel.Coupling(name='CP_Supp_0_xy', controlPoint=RP_Supp_0, surface=N_Supp_0_xy, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=ON)
			myModel.Coupling(name='CP_Supp_0_xyz', controlPoint=RP_Supp_0, surface=N_Supp_0_xyz, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=ON)
			myModel.Coupling(name='CP_Supp_L_x', controlPoint=RP_Supp_L, surface=N_Supp_L_x, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=OFF, u3=OFF, ur1=OFF, ur2=OFF, ur3=ON)
			myModel.Coupling(name='CP_Supp_L_xy', controlPoint=RP_Supp_L, surface=N_Supp_L_xy, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=ON)
			myModel.Coupling(name='CP_Supp_L_xyz', controlPoint=RP_Supp_L, surface=N_Supp_L_xyz, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=ON)
			
		# Torsional actions
		elif Step_n_Mx != 0:
			myModel.Coupling(name='CP_Supp_0_yz', controlPoint=RP_Supp_0, surface=N_Supp_0_yz, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)
			myModel.Coupling(name='CP_Supp_L_yz', controlPoint=RP_Supp_L, surface=N_Supp_L_yz, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)
		
		else:
			myModel.Coupling(name='CP_Supp_0_yz', controlPoint=RP_Supp_0, surface=N_Supp_0_yz, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)
			myModel.Coupling(name='CP_Supp_0_z', controlPoint=RP_Supp_0, surface=N_Supp_0_z, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=OFF, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)
			myModel.Coupling(name='CP_Supp_L_yz', controlPoint=RP_Supp_L, surface=N_Supp_L_yz, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)
			myModel.Coupling(name='CP_Supp_L_z', controlPoint=RP_Supp_L, surface=N_Supp_L_z, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=OFF, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)
	
	# Real supports at lower chord
	else:
		myModel.Coupling(name='CP_Supp_0_xyz', controlPoint=RP_Supp_0, surface=N_Supp_0_xyz, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)
		myModel.Coupling(name='CP_Supp_0_xz', controlPoint=RP_Supp_0, surface=N_Supp_0_xz, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=OFF, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)
		myModel.Coupling(name='CP_Supp_0_z', controlPoint=RP_Supp_0, surface=N_Supp_0_z, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=OFF, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)
		myModel.Coupling(name='CP_Supp_0_y', controlPoint=RP_ForkSupp_0, surface=N_Supp_0_y, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)
		myModel.Coupling(name='CP_Supp_L_xyz', controlPoint=RP_Supp_L, surface=N_Supp_L_xyz, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)
		myModel.Coupling(name='CP_Supp_L_xz', controlPoint=RP_Supp_L, surface=N_Supp_L_xz, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=OFF, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)
		myModel.Coupling(name='CP_Supp_L_z', controlPoint=RP_Supp_L, surface=N_Supp_L_z, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=OFF, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)
		myModel.Coupling(name='CP_Supp_L_y', controlPoint=RP_ForkSupp_L, surface=N_Supp_L_y, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)
		
	''' Bracing '''
	if Brace_type == 'springs+beam' or Brace_type == 'NodestoGround': 
		for intA in range(len(pos_brace)):
			myNodes = []
			[(	myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010010+round(nEL*pos_brace[intA],0)*10000000+intH*10+intB*10000),))), 
				myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10010015+round(nEL*pos_brace[intA],0)*10000000+intH_ZK*10+intB*10000),))), 
				myNodes.append(myAssembly.instances['instbeam'].nodes.sequenceFromLabels(labels=(int(10015010+round(nEL*pos_brace[intA],0)*10000000+intH*10+intB*10000),)))) 
				for intH in range(min_nNH_brace-1, max_nNH_brace) for intH_ZK in range(min_nNH_brace-1, max_nNH_brace-1) for intB in range(nEB+1)]		
			myNameNodes = 'N_Brace_'+str(intA)			
			myAssembly.Set(name=myNameNodes, nodes=myNodes)
			regionbeam = (myAssembly.sets[myNameNodes])
			
			myName = 'CP_Brace_'+str(intA)
			if intA == 0:
				myModel.Coupling(name=myName, controlPoint=RP_Brace_0, surface=regionbeam, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)
			elif intA == 1:
				myModel.Coupling(name=myName, controlPoint=RP_Brace_1, surface=regionbeam, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)
			elif intA == 2:
				myModel.Coupling(name=myName, controlPoint=RP_Brace_2, surface=regionbeam, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)
			elif intA == 3:
				myModel.Coupling(name=myName, controlPoint=RP_Brace_3, surface=regionbeam, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)
			elif intA == 4:
				myModel.Coupling(name=myName, controlPoint=RP_Brace_4, surface=regionbeam, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)
			elif intA == 5:
				myModel.Coupling(name=myName, controlPoint=RP_Brace_5, surface=regionbeam, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)
			elif intA == 6:
				myModel.Coupling(name=myName, controlPoint=RP_Brace_6, surface=regionbeam, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)
			elif intA == 7:
				myModel.Coupling(name=myName, controlPoint=RP_Brace_7, surface=regionbeam, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)
			elif intA == 8:
				myModel.Coupling(name=myName, controlPoint=RP_Brace_8, surface=regionbeam, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)
			elif intA == 9:
				myModel.Coupling(name=myName, controlPoint=RP_Brace_9, surface=regionbeam, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)		
	
	''' Actions '''
	# Axial actions
	if Step_n_Nx != 0 and Load_type_Nx == 'w' and Supp_type == 'real':
		myModel.Coupling(name='CP_Nx_0', controlPoint=RP_Nx_0, surface=N_Nx_0, influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, weightingMethod=UNIFORM, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
		myModel.Coupling(name='CP_Nx_L', controlPoint=RP_Nx_L, surface=N_Nx_L, influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, weightingMethod=UNIFORM, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
	
	if Step_n_Nx != 0 and Load_type_Nx == 'Exp_LTB' and Supp_type == 'real':		
		myModel.Coupling(name='CP_Nx_0', controlPoint=RP_Nx_0, surface=N_Nx_0, influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, weightingMethod=UNIFORM, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
		myModel.Coupling(name='CP_Nx_L', controlPoint=RP_Nx_L, surface=N_Nx_L, influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, weightingMethod=UNIFORM, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)			
	
	# Torsional actions
	if Step_n_Mx != 0:
		myModel.Coupling(name='CP_Mx', controlPoint=RP_Mx, surface=N_Mx, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=ON, ur1=ON, ur2=OFF, ur3=OFF)

	# Actions causing My
	if Step_n_My != 0:
		# Point load - 3-point bending
		if Moment == 'Qz_3P': 
			myModel.Coupling(name='CP_Qz_3P', controlPoint=RP_Qz_3P, surface=N_Qz_3P, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=ON, ur1=ON, ur2=OFF, ur3=OFF)
		
		# Point load - 3-point bending - distributed over an area
		if Moment == 'Qz_3P_area': 
			myModel.Coupling(name='CP_Qz_3P_area_yz', controlPoint=RP_Qz_3P_area, surface=N_Qz_3P_area_yz, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, weightingMethod=UNIFORM, localCsys=None, u1=OFF, u2=ON, u3=ON, ur1=ON, ur2=OFF, ur3=OFF)
			myModel.Coupling(name='CP_Qz_3P_area_z', controlPoint=RP_Qz_3P_area, surface=N_Qz_3P_area_z, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, weightingMethod=UNIFORM, localCsys=None, u1=OFF, u2=OFF, u3=ON, ur1=ON, ur2=OFF, ur3=OFF)
		if Moment == 'Exp_LTB':
			myModel.Coupling(name='CP_Qz_3P_area', controlPoint=RP_Qz_3P_area, surface=N_Qz_3P_area, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, weightingMethod=UNIFORM, localCsys=None, u1=OFF, u2=OFF, u3=ON, ur1=ON, ur2=OFF, ur3=OFF)
			myModel.Coupling(name='CP_Qz_3P_screws', controlPoint=RP_Qz_3P_area, surface=N_Qz_3P_screws, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, weightingMethod=UNIFORM, localCsys=None, u1=OFF, u2=OFF, u3=ON, ur1=ON, ur2=OFF, ur3=OFF)
			myModel.Coupling(name='CP_Qz_3P_area_side', controlPoint=RP_Qz_3P_area, surface=N_Qz_3P_area_side, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, weightingMethod=UNIFORM, localCsys=None, u1=OFF, u2=ON, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)
			
		# Point load - 4-point bending
		if Moment == 'Qz_4P': 
			myModel.Coupling(name='CP_Qz_4P_0', controlPoint=RP_Qz_4P_0, surface=N_Qz_4P_0, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=ON, ur1=ON, ur2=OFF, ur3=OFF)
			myModel.Coupling(name='CP_Qz_4P_L', controlPoint=RP_Qz_4P_L, surface=N_Qz_4P_L, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=ON, ur1=ON, ur2=OFF, ur3=OFF)
		
		# Point load - 4-point bending - distributed over an area
		if Moment == 'Qz_4P_area' or Moment == 'Qz_4P_area_T': 
			myModel.Coupling(name='CP_Qz_4P_area_0', controlPoint=RP_Qz_4P_0, surface=N_Qz_4P_area_0, influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, weightingMethod=UNIFORM, localCsys=None, u1=OFF, u2=ON, u3=ON, ur1=ON, ur2=OFF, ur3=OFF)
			myModel.Coupling(name='CP_Qz_4P_area_L', controlPoint=RP_Qz_4P_L, surface=N_Qz_4P_area_L, influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, weightingMethod=UNIFORM, localCsys=None, u1=OFF, u2=ON, u3=ON, ur1=ON, ur2=OFF, ur3=OFF)
			if Moment == 'Qz_4P_area_T': 
				myModel.Coupling(name='CP_Qz_4P_area_T0', controlPoint=RP_Qz_4P_area_T, surface=RP_Qz_4P_0, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
				myModel.Coupling(name='CP_Qz_4P_area_TL', controlPoint=RP_Qz_4P_area_T, surface=RP_Qz_4P_L, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
				
		# Bending moments My at supports
		if Moment == 'My_const' and Load_type_My == 'w' and Supp_type == 'real':
			myModel.Coupling(name='CP_My_const_0', controlPoint=RP_My_const_0, surface=N_My_const_0, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
			myModel.Coupling(name='CP_My_const_L', controlPoint=RP_My_const_L, surface=N_My_const_L, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
		
	########################################
	''' Loading '''
	########################################	
	if Cal_type == 'LBA':
		maxStep = 1
	elif Cal_type == 'GMNIA':
		maxStep = max(Step_n_Nx,Step_n_Mx,Step_n_My)	# Maximum step number
	
	# Loop over steps
	for Step in range (maxStep):
		# Set basic parameters
		if Step+1 == 1:
			maxnInc = maxnInc_One
			startInc = startInc_One
			minInc = minInc_One
			maxInc = maxInc_One
			StepNameBefore = 'Initial'
			StepName = 'One'
		elif Step+1 == 2:
			maxnInc = maxnInc_Two
			startInc = startInc_Two
			minInc = minInc_Two
			maxInc = maxInc_Two
			StepNameBefore = StepName
			StepName = 'Two'
		elif Step+1 == 3:
			maxnInc = maxnInc_Three
			startInc = startInc_Three
			minInc = minInc_Three
			maxInc = maxInc_Three	
			StepNameBefore = StepName
			StepName = 'Three'
		
		# Create step
		if Cal_type == 'LBA':
			myModel.BuckleStep(name='LBA', previous='Initial', numEigen=6, vectors=12, maxIterations=30)
			StepName = 'LBA'
		elif Cal_type == 'GMNIA':	
			myModel.StaticStep(name=StepName, previous=StepNameBefore, maxNumInc = maxnInc, initialInc = startInc_One, minInc = minInc, maxInc = maxInc, nlgeom=ON, solutionTechnique=Sol_tech)	
		
		# Axial actions
		if Step_n_Nx == Step+1: 
			# Displacement controlled
			if Load_type_Nx == 'w':
				if Supp_type == 'ideal':
					myModel.boundaryConditions['Supp_0'].setValuesInStep(stepName=StepName, u1=wx, ur2=FREED, ur3=FREED)
					myModel.boundaryConditions['Supp_L'].setValuesInStep(stepName=StepName, u1=-wx, ur2=FREED, ur3=FREED)			
				else:
					myModel.DisplacementBC(name='Nx_0', createStepName=StepName, region=RP_Nx_0, u1=wx, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
					myModel.DisplacementBC(name='Nx_L', createStepName=StepName, region=RP_Nx_L, u1=-wx, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
			
			# Force controlled
			elif Load_type_Nx == 'Exp_LTB':
				myModel.ConcentratedForce(name='Nx_0', createStepName=StepName, region=RP_Nx_0, cf1=Nx, distributionType=UNIFORM, field='', localCsys=None)
				myModel.ConcentratedForce(name='Nx_L', createStepName=StepName, region=RP_Nx_L, cf1=-Nx, distributionType=UNIFORM, field='', localCsys=None)
				
			# Force controlled
			else:		
				if Cal_type == 'LBA':
					Sigma_Nx_0 = 1/float(B*H_0)
					Sigma_Nx_L = 1/float(B*H_L)
				elif Cal_type == 'GMNIA':
					Sigma_Nx_0 = Nx/float(B*H_0)
					Sigma_Nx_L = Nx/float(B*H_L)
					
				myModel.SurfaceTraction(name='Nx_0', createStepName=StepName, region=SF_Nx_0, magnitude=Sigma_Nx_0, directionVector=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), distributionType=UNIFORM, field='', localCsys=None, traction=GENERAL, follower=OFF, resultant=ON)
				myModel.SurfaceTraction(name='Nx_L', createStepName=StepName, region=SF_Nx_L, magnitude=Sigma_Nx_L, directionVector=((0.0, 0.0, 0.0), (-1.0, 0.0, 0.0)), distributionType=UNIFORM, field='', localCsys=None, traction=GENERAL, follower=OFF, resultant=ON)
				
				# Bending stress My due to load eccentricity of Nx in z direction
				# At x = 0
				if ez_Nx_0 != 0:
					if Cal_type == 'LBA':
						Sigma_Nx_My_0 = 1*ez_Nx_0/float(H_0**2*B/6)
					elif Cal_type == 'GMNIA':
						Sigma_Nx_My_0 = Nx*ez_Nx_0/float(H_0**2*B/6)
						
					myModel.ExpressionField(name='F_My_0', localCsys=None, description='', expression=' (Z-z_cen_0)/(H_0/2)')		# Defining function for creating a linear stress distribution over the cross-section height
					myModel.SurfaceTraction(name='My_0', createStepName=StepName, region=SF_Nx_0, magnitude=Sigma_Nx_My_0, directionVector=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), distributionType=FIELD, field='F_My_0', localCsys=None, traction=GENERAL, follower=OFF, resultant=ON)
				
				# At x = L
				if ez_Nx_L != 0:
					if Cal_type == 'LBA':
						Sigma_Nx_My_L = 1*ez_Nx_L/float(H_L**2*B/6)	
					elif Cal_type == 'GMNIA':
						Sigma_Nx_My_L = Nx*ez_Nx_L/float(H_L**2*B/6)
						
					myModel.ExpressionField(name='F_My_L', localCsys=None, description='', expression=' (Z-z_cen_L)/(H_L/2)')		# Defining function for creating a linear stress distribution over the cross-section height
					myModel.SurfaceTraction(name='My_L', createStepName=StepName, region=SF_Nx_L, magnitude=Sigma_Nx_My_L, directionVector=((0.0, 0.0, 0.0), (-1.0, 0.0, 0.0)), distributionType=FIELD, field='F_My_L', localCsys=None, traction=GENERAL, follower=OFF, resultant=ON)
				
				# Bending stress Mz due to load eccentricity of Nx in y direction
				# At x = 0
				if ey_Nx_0 != 0:
					if Cal_type == 'LBA':
						Sigma_Nx_Mz_0 = 1*ey_Nx_0/float(B**2*H_0/6)
					elif Cal_type == 'GMNIA':
						Sigma_Nx_Mz_0 = Nx*ey_Nx_0/float(B**2*H_0/6)
						
					myModel.ExpressionField(name='F_Mz_0', localCsys=None, description='', expression=' Y/(B/2)')		# Defining function for creating a linear stress distribution over the cross-section height
					myModel.SurfaceTraction(name='Mz_0', createStepName=StepName, region=SF_Nx_0, magnitude=Sigma_Nx_Mz_0, directionVector=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), distributionType=FIELD, field='F_Mz_0', localCsys=None, traction=GENERAL, follower=OFF, resultant=ON)
				
				# At x = L
				if ey_Nx_L != 0:
					if Cal_type == 'LBA':
						Sigma_Nx_Mz_L = 1*ey_Nx_L/float(B**2*H_L/6)	
					elif Cal_type == 'GMNIA':
						Sigma_Nx_Mz_L = Nx*ey_Nx_L/float(B**2*H_L/6)
						
					myModel.ExpressionField(name='F_Mz_L', localCsys=None, description='', expression=' Y/(B/2)')		# Defining function for creating a linear stress distribution over the cross-section height
					myModel.SurfaceTraction(name='Mz_L', createStepName=StepName, region=SF_Nx_L, magnitude=Sigma_Nx_Mz_L, directionVector=((0.0, 0.0, 0.0), (-1.0, 0.0, 0.0)), distributionType=FIELD, field='F_Mz_L', localCsys=None, traction=GENERAL, follower=OFF, resultant=ON)
			
		# Torsional actions
		if Step_n_Mx == Step+1: 
			# Force controlled
			if Load_type_Mx == 'F':
				if Cal_type == 'LBA':
					myModel.Moment(name='Mx', createStepName=StepName, region=RP_Mx, cm1=1, distributionType=UNIFORM, field='', localCsys=None, follower=ON)
				elif Cal_type == 'GMNIA':
					myModel.Moment(name='Mx', createStepName=StepName, region=RP_Mx, cm1=Mx, distributionType=UNIFORM, field='', localCsys=None, follower=ON)
				
			# Displacement controlled
			else:
				myModel.boundaryConditions['Supp_Mx'].setValuesInStep(stepName=StepName, ur1=rotx)
			
		# Actions causing My
		if Step_n_My == Step+1:
			# Point load - 3-point bending
			if Moment == 'Qz_3P': 
				# Force controlled
				if Load_type_My == 'F':
					if Cal_type == 'LBA':
						myModel.ConcentratedForce(name='Qz_3P', createStepName=StepName, region=RP_Qz_3P, cf3=1, distributionType=UNIFORM, field='', localCsys=None)
					elif Cal_type == 'GMNIA':
						myModel.ConcentratedForce(name='Qz_3P', createStepName=StepName, region=RP_Qz_3P, cf3=F, distributionType=UNIFORM, field='', localCsys=None)
				else:
					myModel.boundaryConditions['Supp_Qz_3P'].setValuesInStep(stepName=StepName, u3=w)
			
			# Point load - 3-point bending - applied as surface traction
			if Moment == 'Qz_3Ps': 
				# Force controlled
				if Load_type_My == 'F':
					if Cal_type == 'LBA':
						Sigma_Qz_3Ps = 1/(B*H_Qz_3Ps)
					elif Cal_type == 'GMNIA':
						Sigma_Qz_3Ps = F/(B*H_Qz_3Ps)
						
					myModel.SurfaceTraction(name='Qz_3Ps', createStepName=StepName, region=SF_Qz_3Ps, magnitude=Sigma_Qz_3Ps, directionVector=((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), distributionType=UNIFORM, field='', localCsys=None, traction=GENERAL, follower=OFF, resultant=ON)
				
			# Point load distributed over an area - 3-point bending
			if Moment == 'Qz_3P_area' or Moment == 'Exp_LTB': 
				# Force controlled
				if Load_type_My == 'F':
					if Cal_type == 'LBA':
						myModel.ConcentratedForce(name='Qz_3P_area', createStepName=StepName, region=RP_Qz_3P_area, cf3=1, distributionType=UNIFORM, field='', localCsys=None)
						if Moment == 'Exp_LTB':
							F_hor = ey_Exp_LTB / ez_Exp_LTB
							myModel.ConcentratedForce(name='F_hor', createStepName=StepName, region=RP_Qz_3P_area, cf2=F_hor, distributionType=UNIFORM, field='', localCsys=None)
					elif Cal_type == 'GMNIA':
						myModel.ConcentratedForce(name='Qz_3P_area', createStepName=StepName, region=RP_Qz_3P_area, cf3=F, distributionType=UNIFORM, field='', localCsys=None)
						if Moment == 'Exp_LTB':
							F_hor = F * ey_Exp_LTB / ez_Exp_LTB
							myModel.ConcentratedForce(name='F_hor', createStepName=StepName, region=RP_Qz_3P_area, cf2=F_hor, distributionType=UNIFORM, field='', localCsys=None)
				else:
					myModel.boundaryConditions['Supp_Qz_3P_area'].setValuesInStep(stepName=StepName, u3=w)
					if Moment == 'Exp_LTB':
						myModel.boundaryConditions['Supp_Qz_3P_area'].setValuesInStep(stepName=StepName, u3=w)
						myModel.TabularAmplitude(name='F_hor_tab', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=data_F_hor_tab)
						myModel.ConcentratedForce(name='F_hor', createStepName=StepName, region=RP_Qz_3P_area, cf2=1.0, amplitude='F_hor_tab', distributionType=UNIFORM, field='', localCsys=None)
					
			# Point load - 4-point bending
			if Moment == 'Qz_4P': 
				# Force controlled
				if Load_type_My == 'F':
					if Cal_type == 'LBA':
						myModel.ConcentratedForce(name='Qz_4P_0', createStepName=StepName, region=RP_Qz_4P_0, cf3=1, distributionType=UNIFORM, field='', localCsys=None)
						myModel.ConcentratedForce(name='Qz_4P_L', createStepName=StepName, region=RP_Qz_4P_L, cf3=1, distributionType=UNIFORM, field='', localCsys=None)
					elif Cal_type == 'GMNIA':
						myModel.ConcentratedForce(name='Qz_4P_0', createStepName=StepName, region=RP_Qz_4P_0, cf3=F, distributionType=UNIFORM, field='', localCsys=None)
						myModel.ConcentratedForce(name='Qz_4P_L', createStepName=StepName, region=RP_Qz_4P_L, cf3=F, distributionType=UNIFORM, field='', localCsys=None)

				# Displacement controlled
				else:
					myModel.boundaryConditions['Supp_Qz_4P_0'].setValuesInStep(stepName=StepName, u3=w)
					myModel.boundaryConditions['Supp_Qz_4P_L'].setValuesInStep(stepName=StepName, u3=w)
			
			# Point load - 4-point bending - applied as surface traction
			if Moment == 'Qz_4Ps': 
				# Force controlled
				if Load_type_My == 'F':
					if Cal_type == 'LBA':
						Sigma_Qz_4Ps_0 = 1/(B*H_Qz_4Ps_0)
						Sigma_Qz_4Ps_L = 1/(B*H_Qz_4Ps_L)
					elif Cal_type == 'GMNIA':
						Sigma_Qz_4Ps_0 = F/(B*H_Qz_4Ps_0)
						Sigma_Qz_4Ps_L = F/(B*H_Qz_4Ps_L)
						
					myModel.SurfaceTraction(name='Qz_4Ps_0', createStepName=StepName, region=SF_Qz_4Ps_0, magnitude=Sigma_Qz_4Ps_0, directionVector=((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), distributionType=UNIFORM, field='', localCsys=None, traction=GENERAL, follower=OFF, resultant=ON)
					myModel.SurfaceTraction(name='Qz_4Ps_L', createStepName=StepName, region=SF_Qz_4Ps_L, magnitude=Sigma_Qz_4Ps_L, directionVector=((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), distributionType=UNIFORM, field='', localCsys=None, traction=GENERAL, follower=OFF, resultant=ON)
			
			# Point load distributed over an area - 4-point bending
			if Moment == 'Qz_4P_area': 
				# Force controlled
				if Load_type_My == 'F':
					if Cal_type == 'LBA':
						myModel.ConcentratedForce(name='Qz_4P_area_0', createStepName=StepName, region=RP_Qz_4P_0, cf3=1, distributionType=UNIFORM, field='', localCsys=None)
						myModel.ConcentratedForce(name='Qz_4P_area_L', createStepName=StepName, region=RP_Qz_4P_L, cf3=1, distributionType=UNIFORM, field='', localCsys=None)
					elif Cal_type == 'GMNIA':
						myModel.ConcentratedForce(name='Qz_4P_area_0', createStepName=StepName, region=RP_Qz_4P_0, cf3=F, distributionType=UNIFORM, field='', localCsys=None)
						myModel.ConcentratedForce(name='Qz_4P_area_L', createStepName=StepName, region=RP_Qz_4P_L, cf3=F, distributionType=UNIFORM, field='', localCsys=None)
				else:
					myModel.boundaryConditions['Supp_Qz_4P_area_0'].setValuesInStep(stepName=StepName, u3=w)
					myModel.boundaryConditions['Supp_Qz_4P_area_L'].setValuesInStep(stepName=StepName, u3=w)
			
			# Point load distributed over an area - 4-point bending with traverse
			if Moment == 'Qz_4P_area_T': 
				# Force controlled
				if Load_type_My == 'F':
					if Cal_type == 'LBA':
						myModel.ConcentratedForce(name='Qz_4P_area_T', createStepName=StepName, region=RP_Qz_4P_area_T, cf3=1, distributionType=UNIFORM, field='', localCsys=None)
					elif Cal_type == 'GMNIA':
						myModel.ConcentratedForce(name='Qz_4P_area_T', createStepName=StepName, region=RP_Qz_4P_area_T, cf3=F, distributionType=UNIFORM, field='', localCsys=None)
				else:
					myModel.boundaryConditions['Supp_Qz_4P_area_T'].setValuesInStep(stepName=StepName, u3=w)
					
			# Constant moment My
			if Moment == 'My_const': 
				# Force controlled
				if Load_type_My == 'F':
					# At x = 0
					if Cal_type == 'LBA':
						Sigma_My_const_0 = 1/float(H_0**2*B/6)
					elif Cal_type == 'GMNIA':
						Sigma_My_const_0 = F/float(H_0**2*B/6)
						
					myModel.ExpressionField(name='F_My_const_0', localCsys=None, description='', expression=' (Z-z_cen_0)/(H_0/2)')		# Defining function for creating a linear stress distribution over the cross-section height
					myModel.SurfaceTraction(name='My_const_0', createStepName=StepName, region=SF_My_const_0, magnitude=Sigma_My_const_0, directionVector=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), distributionType=FIELD, field='F_My_const_0', localCsys=None, traction=GENERAL, follower=OFF, resultant=ON)
					
					# At x = L
					if Cal_type == 'LBA':
						Sigma_My_const_L = 1/float(H_L**2*B/6)	
					elif Cal_type == 'GMNIA':
						Sigma_My_const_L = F/float(H_L**2*B/6)
						
					myModel.ExpressionField(name='F_My_const_L', localCsys=None, description='', expression=' (Z-z_cen_L)/(H_L/2)')		# Defining function for creating a linear stress distribution over the cross-section height
					myModel.SurfaceTraction(name='My_const_L', createStepName=StepName, region=SF_My_const_L, magnitude=Sigma_My_const_L, directionVector=((0.0, 0.0, 0.0), (-1.0, 0.0, 0.0)), distributionType=FIELD, field='F_My_const_L', localCsys=None, traction=GENERAL, follower=OFF, resultant=ON)
					
				# Displacement controlled	
				else:
					if Supp_type == 'ideal':
						myModel.boundaryConditions['Supp_0'].setValuesInStep(stepName=StepName, ur2=w)
						myModel.boundaryConditions['Supp_L'].setValuesInStep(stepName=StepName, ur2=-w)
					else:
						myModel.DisplacementBC(name='My_const_0', createStepName=StepName, region=RP_My_const_0, u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=w, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
						myModel.DisplacementBC(name='My_const_L', createStepName=StepName, region=RP_My_const_L, u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=-w, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
						
			# Constant line load qz
			if Moment == 'qz': 
				# Force controlled
				if Load_type_My == 'F':
					if Cal_type == 'LBA':
						Sigma_qz = 1/B	
					elif Cal_type == 'GMNIA':
						Sigma_qz = F/B
						
					myModel.SurfaceTraction(name='qz', createStepName=StepName, region=SF_qz, magnitude=Sigma_qz, directionVector=((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), distributionType=UNIFORM, field='', localCsys=None, traction=GENERAL, follower=OFF, resultant=ON)
			
	########################################
	''' Output settings '''
	########################################
	# Add results of further Steps
	if Cal_type == 'GMNIA':
		myModel.FieldOutputRequest(name='F-Output-1', createStepName='One', variables=('CDISP', 'CF', 'CSTRESS', 'EEQUT', 'LE', 'MVF', 'NFORCSO', 'PE', 'PEEQ', 'PEMAG',  'RF', 'S', 'SDV', 'SEQUT', 'TE', 'TEEQ', 'TEVOL', 'U'))
		myModel.HistoryOutputRequest(name='H-Output-1', createStepName='One', variables=PRESELECT)
		
	if maxStep >= 2: 
		myModel.FieldOutputRequest(name='F-Output-2', createStepName='Two', variables=('CDISP', 'CF', 'CSTRESS', 'EEQUT', 'LE', 'MVF', 'NFORCSO', 'PE', 'PEEQ', 'PEMAG',  'RF', 'S', 'SDV', 'SEQUT', 'TE', 'TEEQ', 'TEVOL', 'U'))
		myModel.HistoryOutputRequest(name='H-Output-2', createStepName='Two', variables=PRESELECT)
	if maxStep >= 3: 
		myModel.FieldOutputRequest(name='F-Output-3', createStepName='Three', variables=('CDISP', 'CF', 'CSTRESS', 'EEQUT', 'LE', 'MVF', 'NFORCSO', 'PE', 'PEEQ', 'PEMAG',  'RF', 'S',  'SDV', 'SEQUT', 'TE', 'TEEQ', 'TEVOL', 'U'))
		myModel.HistoryOutputRequest(name='H-Output-3', createStepName='Three', variables=PRESELECT)
		
	# Add results of springs
	if Cal_type == 'GMNIA':
		if Step_n_Nx != 0 and M_mu_supp != 0:
			myModel.FieldOutputRequest(name='Connectors-1', createStepName='One', variables=('CEF', 'CUE'), region=myAssembly.sets['Wire_Spring_0'], sectionPoints=DEFAULT, rebar=EXCLUDE)
			if maxStep >= 2: myModel.FieldOutputRequest(name='Connectors-2', createStepName='Two', variables=('CEF', 'CUE'), region=myAssembly.sets['Wire_Spring_0'], sectionPoints=DEFAULT, rebar=EXCLUDE)
			if maxStep >= 3: myModel.FieldOutputRequest(name='Connectors-2', createStepName='Three', variables=('CEF', 'CUE'), region=myAssembly.sets['Wire_Spring_0'], sectionPoints=DEFAULT, rebar=EXCLUDE)
				
		if Brace_type != 'none':
			for intA in range(len(pos_brace)):
				myName = tuple(['SP_Brace_'+str(intA)])
				myNameOutput = 'H-Output-'+str(intA+1)
				mdb.models[modelName].HistoryOutputRequest(name=myNameOutput, createStepName='One', variables=('S11', 'E11'), springs=(myName), sectionPoints=DEFAULT, rebar=EXCLUDE)
				if maxStep >= 2: mdb.models[modelName].HistoryOutputRequest(name=myNameOutput, createStepName='Two', variables=('S11', 'E11'), springs=(myName), sectionPoints=DEFAULT, rebar=EXCLUDE)
				if maxStep >= 3: mdb.models[modelName].HistoryOutputRequest(name=myNameOutput, createStepName='Three', variables=('S11', 'E11'), springs=(myName), sectionPoints=DEFAULT, rebar=EXCLUDE)
			
	# Add internal forces of equivalent beam bracing
	if Brace_type == 'springs+beam':
		myModel.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'SDV', 'SEQUT', 'PE', 'PEEQ', 'PEMAG', 'LE', 'TE', 'TEEQ', 'TEVOL', 'EEQUT', 'U', 'RF', 'CF', 'NFORCSO', 'CSTRESS', 'CDISP', 'MVF'))

	# Settings View
	session.viewports['Viewport: 1'].view.rotate(xAngle=135, yAngle=0, zAngle=45, mode=TOTAL)

	################################################################################################################################################################################################################################################
	''' Create & Submit Job '''
	################################################################################################################################################################################################################################################
	# Create Job
	if Cal_type == 'GMNIA':
		jobName = 'Cal_'+str(i)+'_'+Cal_type
	elif Cal_type == 'LBA':
		jobName = 'Cal_'+str(i)+'_'+Cal_type
		
	myJob = mdb.Job(name=jobName, model=myModel, type=ANALYSIS, memoryUnits=PERCENTAGE, memory=90, numCpus=4, numDomains=4)
	if Material_law == 'UMAT' and UMAT_format == 'for':
		myJob.setValues(userSubroutine='UMAT.for')
	elif Material_law == 'UMAT' and UMAT_format == 'f':
		myJob.setValues(userSubroutine='UMAT.f')
	
	if submit_job == 1:
		# Submit Job / Start computation
		try:
			myJob.submit(datacheckJob=False)
			myJob.waitForCompletion()
			conv = True
		except AbaqusException, message:
			conv = False
		
		################################################################################################################################################################################################################################################
		''' Postprocessing '''
		################################################################################################################################################################################################################################################
		''' Output '''
		####################################	
		myOdb = str(jobName)+'.odb'
		session.viewports['Viewport: 1'].setValues(displayedObject=session.openOdb(name=myOdb))

		# Open the output database
		odb = openOdb(myOdb)

		# Name of .txt with results
		saveName = '_Output'

		# Add Output here
			
		# Close odb
		odb.close()
		
	####################################
	''' Fin '''
	####################################
	# Delete default model
	if i == 0:
		mdb.Model(modelType=STANDARD_EXPLICIT, name='Model-1')	# Create default model if not already existing
		del mdb.models['Model-1']								# Delete default model
	
	# Safe model
	if safe_model == 1: mdb.saveAs(pathName=modelName_save)
	
	# Delete job
	if del_job == 1:
		del mdb.jobs[jobName]
		
	# Clear calculation files in work_directory
	if del_files == 1:
		myfile = [jobName+".com", jobName+".dat", jobName+".inp", jobName+".ipm", jobName+".log", jobName+".msg", jobName+".prt", jobName+".sta", jobName+".odb"]
		for intF in range(len(myfile)):
			if os.path.isfile(myfile[intF]):
				try:
					os.remove(myfile[intF])	
				except OSError:
					pass
		files = os.listdir('.')				# List all files in the current directory
		for file in files:					# Move files
			if file.endswith('.pyc'):
				try:
					os.remove(file)	
				except OSError:
					pass
################################################################################################################################################################################################################################################
''' Fin/End/Ende '''
################################################################################################################################################################################################################################################
