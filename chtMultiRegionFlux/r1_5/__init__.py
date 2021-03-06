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
def main_standalone( argc, argv ):

    from Foam.OpenFOAM.include import setRootCase
    args = setRootCase( argc, argv )

    from Foam.OpenFOAM.include import createTime
    runTime = createTime( args )    
    
    from regionProperties import regionProperties
    rp = regionProperties( runTime )

    from fluid import createFluidMeshes
    fluidRegions = createFluidMeshes( rp, runTime )
    
    from solid import createSolidMeshes
    solidRegions = createSolidMeshes( rp, runTime )
    
    from fluid import createFluidFields
    pdf, thermof, rhof, Kf, Uf, phif, turb, DpDtf, ghf, initialMassf, pRef = createFluidFields( fluidRegions, runTime, rp )
    
    from solid import createSolidField
    rhos, cps, rhosCps, Ks, Ts = createSolidField( solidRegions, runTime )
    
    from Foam.finiteVolume.cfdTools.general.include import initContinuityErrs
    cumulativeContErr = initContinuityErrs()

    from Foam.finiteVolume.cfdTools.general.include import readTimeControls
    adjustTimeStep, maxCo, maxDeltaT = readTimeControls( runTime )
    
    if fluidRegions.size() :
       from fluid import compressubibleMultiRegionCourantNo
       CoNum = compressubibleMultiRegionCourantNo( fluidRegions, runTime, rhof, phif )
                
       from fluid import setInitialDeltaT
       runTime = setInitialDeltaT( runTime, adjustTimeStep, maxCo, maxDeltaT, CoNum )
       pass
    
    from Foam.OpenFOAM import ext_Info, nl
    
    while runTime.run():
       from Foam.finiteVolume.cfdTools.general.include import readTimeControls
       adjustTimeStep, maxCo, maxDeltaT = readTimeControls( runTime )
       
       if fluidRegions.size() :
          from fluid import compressubibleMultiRegionCourantNo
          CoNum = compressubibleMultiRegionCourantNo( fluidRegions, runTime, rhof, phif )

          from Foam.finiteVolume.cfdTools.general.include import setDeltaT   
          runTime = setDeltaT( runTime, adjustTimeStep, maxCo, maxDeltaT, CoNum )
          pass
       
       runTime.increment()
       ext_Info()<< "Time = " << runTime.timeName() << nl << nl
       
       for i in range( fluidRegions.size() ):
           ext_Info() << "\nSolving for fluid region " << fluidRegions[ i ].name() << nl
           
           from fluid import readFluidMultiRegionPISOControls
           piso, nCorr, nNonOrthCorr, momentumPredictor, transonic, nOuterCorr = readFluidMultiRegionPISOControls( fluidRegions[ i ] )

           from fluid import solveFluid
           cumulativeContErr = solveFluid( i, fluidRegions, pdf, thermof, rhof, Kf, Uf, phif, turb, DpDtf, ghf, initialMassf, pRef,\
                                           nCorr, nNonOrthCorr, momentumPredictor, transonic, nOuterCorr, cumulativeContErr )
           
           pass
        
       for i in range( solidRegions.size() ):
           ext_Info() << "\nSolving for solid region " << solidRegions[ i ].name() << nl
           
           from solid import readSolidMultiRegionPISOControls
           piso, nNonOrthCorr = readSolidMultiRegionPISOControls( solidRegions[ i ] )
               
           from solid import solveSolid
           solveSolid( i, rhosCps,  Ks, Ts, nNonOrthCorr )
           pass
       
       runTime.write();

       ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" \
            << "  ClockTime = " << runTime.elapsedClockTime() << " s" \
            << nl << nl
       

    ext_Info() << "End\n"
    pass

       
    import os
    return os.EX_OK

    
#--------------------------------------------------------------------------------------
argv = None
import sys, os

from Foam import FOAM_VERSION
if FOAM_VERSION( "==", "010500" ):
    if __name__ == "__main__" :
        argv = sys.argv
        os._exit( main_standalone( len( argv ), argv ) )
        pass
else:
    from Foam.OpenFOAM import ext_Info
    ext_Info() << "\n\n To use this solver, it is necessary to SWIG OpenFOAM-1.5 \n"    
    pass

