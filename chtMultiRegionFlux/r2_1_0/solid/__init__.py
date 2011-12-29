#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV
##


#-----------------------------------------------------------------
from Foam import ref, man


#-----------------------------------------------------------------
def createSolidMeshes( rp, runTime ):
    
    solidRegions = list()
    for index in range( rp.solidRegionNames().size() ):
        ref.ext_Info()<< "Create solid mesh for region " << rp.solidRegionNames()[ index ] \
            << " for time = " << runTime.timeName() << ref.nl << ref.nl
        
        solidRegions.append( man.fvMesh( man.IOobject ( rp.solidRegionNames()[ index ],
                                                        ref.fileName( runTime.timeName() ),
                                                        runTime,
                                                        ref.IOobject.MUST_READ ) ) )
        pass

    return solidRegions


#---------------------------------------------------------------------
def createSolidField( solidRegions, runTime ):

    thermos = list()
    for index in range( solidRegions.__len__() ):
        ref.ext_Info() << "*** Reading solid mesh thermophysical properties for region " \
                       << solidRegions[ index ].name() << ref.nl << ref.nl

        ref.ext_Info() << "    Adding to thermos\n" << ref.nl
        thermos.append( man.basicSolidThermo.New( solidRegions[ index ] ) )
        pass
   
    return thermos


#-----------------------------------------------------------------------------------------------------
def readSolidMultiRegionPIMPLEControls( mesh ):
    pimple = mesh.solutionDict().subDict( ref.word( "PIMPLE" ) )
    nNonOrthCorr = pimple.lookupOrDefault( ref.word( "nNonOrthogonalCorrectors" ), 0 )
    
    return pimple, nNonOrthCorr


#-------------------------------------------------------------------------------------------------------
def readSolidTimeControls( runTime ):
    maxDi = runTime.controlDict().lookupOrDefault( ref.word( "maxDi" ) , 10.0 )
    
    return maxDi


#-------------------------------------------------------------------------------------------------------
def solidRegionDiffNo( mesh, runTime, Cprho, kappa ):
    DiNum = 0.0
    meanDiNum = 0.0

    #- Take care: can have fluid domains with 0 cells so do not test for
    #  zero internal faces.
    kapparhoCpbyDelta = mesh.deltaCoeffs() * ref.fvc.interpolate( kappa ) / ref.fvc.interpolate(Cprho)
    DiNum = kapparhoCpbyDelta.internalField().gMax() * runTime.deltaT().value()
    meanDiNum = kapparhoCpbyDelta.average().value() * runTime.deltaT().value()
    
    ref.ext_Info() << "Region: " << mesh.name() << " Diffusion Number mean: " << meanDiNum << " max: " << DiNum << ref.nl

    return DiNum


#-------------------------------------------------------------------------------------------------------
def solidRegionDiffusionNo( solidRegions, runTime, thermos ):
    DiNum = -ref.GREAT
    
    for index in range( solidRegions.__len__() ):
        mesh, thermo, rho, cp, kappa, T = setRegionSolidFields( index, solidRegions, thermos )
        DiNum = max( solidRegionDiffNo( solidRegions[ index ], runTime, rho * cp, kappa ), DiNum )
        pass
    
    return DiNum

#-------------------------------------------------------------------------------------------------------
def setRegionSolidFields( i, solidRegions, thermos ):
    mesh = solidRegions[ i ]
    thermo = thermos[ i ]

    rho = man.volScalarField( thermo.rho(), man.Deps( thermo ) )

    cp = man.volScalarField( thermo.Cp(), man.Deps( thermo ) )
    
    kappa = man.volScalarField( thermo.K(), man.Deps( thermo ) )
    # tmp<volSymmTensorField> tK = thermo.directionalK();
    
    # const volSymmTensorField& K = tK();
    T = man.volScalarField( thermo.T(), man.Deps( thermo ) )
    
    return mesh, thermo, rho, cp, kappa, T


#-------------------------------------------------------------------------------------------------------
def solveSolid( mesh, thermo, rho, cp, kappa, T, nNonOrthCorr, finalIter ):
    if finalIter:
        mesh.add( ref.word( "finalIteration" ), True )
        pass

    for index in range( nNonOrthCorr + 1 ):
       TEqn = ref.fvm.ddt( rho * cp, T ) - ref.fvm.laplacian( kappa, T )
       TEqn.relax()
       TEqn.solve( mesh.solver( T.select( finalIter ) ) )
       pass
  
    ref.ext_Info()<< "Min/max T:" << T.ext_min() << ' ' << T.ext_max() << ref.nl
    
    thermo.correct()

    if finalIter:
        mesh.remove( ref.word( "finalIteration" ) )
        pass
    pass
