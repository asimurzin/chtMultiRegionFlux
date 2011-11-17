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


#------------------------------------------------------------------------------------
from Foam import ref, man

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def readPIMPLEControls( runTime ):

    solutionDict = ref.fvSolution( runTime )
    
    pimple = solutionDict.subDict( ref.word( "PIMPLE" ) )
    nOuterCorr = pimple.lookupOrDefault( ref.word( "nOuterCorrectors" ), 1);
    
    return nOuterCorr


#--------------------------------------------------------------------------------------
def setInitialMultiRegionDeltaT( adjustTimeStep, runTime, CoNum, DiNum, maxCo, maxDi, maxDeltaT ):
    if adjustTimeStep:
       if runTime.timeIndex() == 0 and ( CoNum > ref.SMALL or DiNum > ref.SMALL ) :
          if CoNum < ref.SMALL:
             CoNum = ref.SMALL
             pass
          if DiNum < ref.SMALL:
             DiNum = ref.SMALL
             pass
          runTime.setDeltaT( min( min( maxCo / CoNum, maxDi / DiNum ) * runTime.deltaT().value(), maxDeltaT ) )
          ref.ext_Info() << "deltaT = " <<  runTime.deltaT().value() << ref.nl
          pass
    return runTime, CoNum, DiNum


#--------------------------------------------------------------------------------------
def setMultiRegionDeltaT( adjustTimeStep, runTime, CoNum, DiNum, maxCo, maxDi, maxDeltaT ):
    if adjustTimeStep:
       if CoNum == -ref.GREAT:
          CoNum = ref.SMALL
          pass
       if DiNum == -ref.GREAT:
          DiNum = ref.SMALL
          pass
       
       maxDeltaTFluid = maxCo / ( CoNum + ref.SMALL )
       maxDeltaTSolid = maxDi / ( DiNum + ref.SMALL )
       
       deltaTFluid = min( min( maxDeltaTFluid, 1.0 + 0.1 * maxDeltaTFluid ), 1.2 );
       runTime.setDeltaT( min( min( deltaTFluid, maxDeltaTSolid ) * runTime.deltaT().value(), maxDeltaT ) )
       ref.ext_Info() << "deltaT = " <<  runTime.deltaT().value() << ref.nl
       pass
    return runTime, CoNum, DiNum


#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )
    
    rp = ref.compressible.regionProperties( runTime )
    
    from fluid import createFluidMeshes
    fluidRegions = createFluidMeshes( rp, runTime )
    
    from solid import createSolidMeshes,createSolidField
    solidRegions = createSolidMeshes( rp,runTime )
    

    from fluid import createFluidFields
    thermoFluid, rhoFluid, KFluid, UFluid, phiFluid, gFluid, turbulence, DpDtFluid, \
                           initialMassFluid, ghFluid, ghfFluid, p_rghFluid, radiation = createFluidFields( fluidRegions, runTime )

    from solid import createSolidField
    thermos = createSolidField( solidRegions, runTime )
    
    from fluid import initContinuityErrs
    cumulativeContErr = initContinuityErrs( fluidRegions )
    
    adjustTimeStep, maxCo, maxDeltaT = ref.readTimeControls( runTime )
    
    from solid import readSolidTimeControls
    maxDi= readSolidTimeControls( runTime )
    
    from fluid import compressubibleMultiRegionCourantNo
    CoNum = compressubibleMultiRegionCourantNo( fluidRegions, runTime, rhoFluid, phiFluid )
    
    from solid import solidRegionDiffusionNo
    DiNum = solidRegionDiffusionNo( solidRegions, runTime, thermos )
    
    runTime, CoNum, DiNum = setInitialMultiRegionDeltaT( adjustTimeStep, runTime, CoNum, DiNum, maxCo, maxDi, maxDeltaT )

    while runTime.run() :
        adjustTimeStep, maxCo, maxDeltaT = ref.readTimeControls(runTime)
        
        maxDi= readSolidTimeControls( runTime )
        
        nOuterCorr = readPIMPLEControls( runTime )
        
        CoNum = compressubibleMultiRegionCourantNo( fluidRegions, runTime, rhoFluid, phiFluid )

        DiNum = solidRegionDiffusionNo( solidRegions, runTime, thermos )
           
        runTime, CoNum, DiNum = setMultiRegionDeltaT( adjustTimeStep, runTime, CoNum, DiNum, maxCo, maxDi, maxDeltaT )
                
        runTime.increment()
        
        ref.ext_Info()<< "Time = " << runTime.timeName() << ref.nl << ref.nl
                
        if nOuterCorr != 1 :
            for i in range( fluidRegions.__len__() ):
                from fluid import setRegionFluidFields
                mesh, thermo, rho, K, U, phi, turb, DpDt, p, psi, h, initialMass, p_rgh, gh, ghf, rad = \
                                     setRegionFluidFields( i, fluidRegions, thermoFluid, rhoFluid, KFluid, UFluid, \
                                                           phiFluid, turbulence, DpDtFluid, initialMassFluid, ghFluid, ghfFluid, p_rghFluid, radiation )
                
                from fluid import storeOldFluidFields
                storeOldFluidFields( p, rho )
                pass
            pass
        
        # --- PIMPLE loop
        for oCorr in range( nOuterCorr ):
            finalIter = ( oCorr == nOuterCorr-1 )
            for i in range( fluidRegions.__len__() ):
                ref.ext_Info() << "\nSolving for fluid region " << fluidRegions[ i ].name() << ref.nl

                from fluid import setRegionFluidFields
                mesh, thermo, rho, K, U, phi, turb, DpDt, p, psi, h, initialMass, p_rgh, gh, ghf, rad = \
                      setRegionFluidFields( i, fluidRegions, thermoFluid, rhoFluid, KFluid, UFluid, \
                                            phiFluid, turbulence, DpDtFluid, initialMassFluid, ghFluid, ghfFluid, p_rghFluid, radiation )
                
                from fluid import readFluidMultiRegionPIMPLEControls
                pimple, nCorr, nNonOrthCorr, momentumPredictor = readFluidMultiRegionPIMPLEControls( mesh ) 
                
                from fluid import solveFluid
                cumulativeContErr = solveFluid( i, mesh, thermo, rad, thermoFluid, rho, K, U, phi, h, turb, DpDt, p, psi, initialMass, p_rgh, gh, \
                                                ghf, oCorr, nCorr, nOuterCorr, nNonOrthCorr, momentumPredictor, cumulativeContErr, finalIter )
                
                pass
                
            for i in range( solidRegions.__len__() ):
               ref.ext_Info() << "\nSolving for solid region " << solidRegions[ i ].name() << ref.nl
               
               from solid import setRegionSolidFields
               mesh, thermo, rho, cp, K, T = setRegionSolidFields( i, solidRegions, thermos )
               
               from solid import readSolidMultiRegionPIMPLEControls
               pimple, nNonOrthCorr = readSolidMultiRegionPIMPLEControls( mesh )
               
               from solid import solveSolid
               solveSolid( mesh, thermo, rho, cp, K, T, nNonOrthCorr, finalIter )
               pass                
            pass
        pass
        runTime.write()

        ref.ext_Info()<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" \
            << "  ClockTime = " << runTime.elapsedClockTime() << " s" \
            << ref.nl << ref.nl    

    ref.ext_Info() << "End\n"
    
    import os
    return os.EX_OK

    
#--------------------------------------------------------------------------------------
argv = None
import sys, os
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020000" ):
    if __name__ == "__main__" :
        argv = sys.argv
        os._exit( main_standalone( len( argv ), argv ) )
        pass
else:
    ref.ext_Info() << "\n\n To use this solver, it is necessary to SWIG OpenFOAM-2.0.0 or higher \n"    
    pass


#--------------------------------------------------------------------------------------

