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
def createFluidMeshes( rp, runTime ) :
    
    fluidRegions = list()
    for index in range( rp.fluidRegionNames().size() ) :
        ref.ext_Info()<< "Create fluid mesh for region " << rp.fluidRegionNames()[ index ] \
                  << " for time = " << runTime.timeName() << ref.nl << ref.nl
        mesh = man.fvMesh( man.IOobject( rp.fluidRegionNames()[ index ],
                                         ref.fileName( runTime.timeName() ),
                                         runTime,
                                         ref.IOobject.MUST_READ ) )
        fluidRegions.append( mesh )
        pass
    
    return fluidRegions
    
    
#-------------------------------------------------------------------
def createFluidFields( fluidRegions, runTime ) :
    
    # Initialise fluid field pointer lists
    thermoFluid = list() 
    rhoFluid = list()
    kappaFluid = list()
    UFluid = list()
    phiFluid = list()
    gFluid = list()
    turbulence =list()
    p_rghFluid = list()
    ghFluid = list()
    ghfFluid = list()
    radiation =list()
    KFluid = list()
    dpdtFluid =list()
    initialMassFluid = list()
    
    #Populate fluid field pointer lists

    for index in range( fluidRegions.__len__() ) :
        ref.ext_Info() << "*** Reading fluid mesh thermophysical properties for region " \
            << fluidRegions[ index ].name() << ref.nl << ref.nl

        ref.ext_Info()<< "    Adding to thermoFluid\n" << ref.nl
        
        thermo = man.basicRhoThermo.New( fluidRegions[ index ] )
        thermoFluid.append( thermo )
        
        ref.ext_Info()<< "    Adding to rhoFluid\n" << ref.nl
        rhoFluid.append( man.volScalarField( man.IOobject( ref.word( "rho" ), 
                                                           ref.fileName( runTime.timeName() ), 
                                                           fluidRegions[ index ], 
                                                           ref.IOobject.NO_READ, 
                                                           ref.IOobject.AUTO_WRITE ),
                                              man.volScalarField( thermoFluid[ index ].rho(), man.Deps( thermoFluid[ index ] ) ) ) )
        
        ref.ext_Info()<< "    Adding to kappaFluid\n" << ref.nl
        kappaFluid.append( man.volScalarField( man.IOobject( ref.word( "kappa" ),
                                                             ref.fileName( runTime.timeName() ),
                                                             fluidRegions[ index ],
                                                             ref.IOobject.NO_READ,
                                                             ref.IOobject.NO_WRITE ),
                                               man.volScalarField( thermoFluid[ index ].Cp() * thermoFluid[ index ].alpha(), 
                                                                   man.Deps( thermoFluid[ index ] ) ) ) )
                                                       
        ref.ext_Info()<< "    Adding to UFluid\n" << ref.nl
        UFluid.append( man.volVectorField( man.IOobject( ref.word( "U" ),
                                                         ref.fileName( runTime.timeName() ),
                                                         fluidRegions[ index ],
                                                         ref.IOobject.MUST_READ,
                                                         ref.IOobject.AUTO_WRITE ),
                                           fluidRegions[ index ] ) )
        
        ref.ext_Info()<< "    Adding to phiFluid\n" << ref.nl
        phiFluid.append( man.surfaceScalarField( man.IOobject( ref.word( "phi" ),
                                                               ref.fileName( runTime.timeName() ),
                                                               fluidRegions[ index ],
                                                               ref.IOobject.READ_IF_PRESENT,
                                                               ref.IOobject.AUTO_WRITE),
                                                  man.linearInterpolate( rhoFluid[ index ] * UFluid[ index ] ) & 
                                                  man.surfaceVectorField( fluidRegions[ index ].Sf(), man.Deps( fluidRegions[ index ] ) ) ) )
        
        ref.ext_Info()<< "    Adding to gFluid\n" << ref.nl
        gFluid.append( man.uniformDimensionedVectorField( man.IOobject( ref.word( "g" ),
                                                                        ref.fileName( runTime.constant() ),
                                                                        fluidRegions[ index ],
                                                                        ref.IOobject.MUST_READ,
                                                                        ref.IOobject.NO_WRITE ) ) )        
        
        ref.ext_Info()<< "    Adding to turbulence\n" << ref.nl
        turbulence.append( man.compressible.turbulenceModel.New( rhoFluid[ index ],
                                                                 UFluid[ index ],
                                                                 phiFluid[ index ],
                                                                 thermoFluid[ index ] ) )
        ref.ext_Info() << "    Adding to ghFluid\n" << ref.nl
        ghFluid.append( man.volScalarField( ref.word( "gh" ) , 
                                            gFluid[ index ] & man.volVectorField( fluidRegions[ index ].C(), man.Deps( fluidRegions[ index ] ) ) ) )

        ref.ext_Info() << "    Adding to ghfFluid\n" << ref.nl
        ghfFluid.append( man.surfaceScalarField( ref.word( "ghf" ), 
                                                 gFluid[ index ] & man.surfaceVectorField( fluidRegions[ index ].Cf(), man.Deps( fluidRegions[ index ] ) ) ) )

        p_rghFluid.append( man.volScalarField( man.IOobject( ref.word( "p_rgh" ),
                                                             ref.fileName( runTime.timeName() ),
                                                             fluidRegions[ index ],
                                                             ref.IOobject.MUST_READ,
                                                             ref.IOobject.AUTO_WRITE ),
                                               fluidRegions[ index ] ) )
        # Force p_rgh to be consistent with p
        p_rghFluid[ index ] << thermoFluid[ index ].p()() - rhoFluid[ index ] * ghFluid[ index ]
        
        radiation.append( man.radiation.radiationModel.New( man.volScalarField( thermoFluid[ index ].T(), man.Deps( thermoFluid[ index ] ) ) ) )
        
        initialMassFluid.append( ref.fvc.domainIntegrate( rhoFluid[ index ] ).value()  )
        
        ref.ext_Info() << "    Adding to KFluid\n" << ref.nl
        KFluid.append( man.volScalarField( ref.word( "K" ),
                                           man.volScalarField( 0.5 * UFluid[ index ].magSqr(), man.Deps( UFluid[ index ] ) ) ) )
        
        ref.ext_Info()<< "    Adding to dpdtFluid\n" << ref.nl
        dpdtFluid.append( man.volScalarField( ref.word( "dpdt" ), 
                                              man.volScalarField( ref.fvc.ddt( thermoFluid[ index ].p() ), man.Deps( thermoFluid[ index ] ) ) ) )

    
    return thermoFluid, rhoFluid, kappaFluid, UFluid, phiFluid, gFluid, turbulence, KFluid, dpdtFluid, initialMassFluid, ghFluid, ghfFluid, p_rghFluid, radiation
        

#-----------------------------------------------------------------------------------------------------------------------
def compressibleCourantNo( mesh, runTime, rho, phi ) :
    
    CoNum = 0.0
    meanCoNum = 0.0
    tmp = ref.fvc.surfaceSum( phi.mag() )
    sumPhi = tmp.internalField() / rho.internalField() 

    CoNum = 0.5 * ( sumPhi / mesh.V().field() ).gMax() * runTime.deltaTValue()

    meanCoNum = 0.5 * ( sumPhi.gSum() / mesh.V().field().gSum() ) * runTime.deltaTValue()

    ref.ext_Info() << "Region: " << mesh.name() << " Courant Number mean: " << meanCoNum \
        << " max: " << CoNum << ref.nl
    
    return  CoNum
            

#-------------------------------------------------------------------------------------------------------------------------
def compressubibleMultiRegionCourantNo(fluidRegions, runTime, rhoFluid, phiFluid ) :
    
    CoNum = -ref.GREAT
    for index in range( fluidRegions.__len__() ):
        CoNum = max( CoNum, compressibleCourantNo( fluidRegions[ index ], runTime, rhoFluid[ index ], phiFluid[ index ] ) )
        pass
    
    return CoNum
    
    
#-------------------------------------------------------------------------------------------------------------------------
def readFluidMultiRegionPIMPLEControls( mesh ) :
    
    pimple = mesh.solutionDict().subDict( ref.word( "PIMPLE" ) )
    nCorr = pimple.lookupOrDefault( ref.word( "nCorrectors" ), 1);
    nNonOrthCorr = pimple.lookupOrDefault( ref.word( "nNonOrthogonalCorrectors" ), 0 )
    momentumPredictor =pimple.lookupOrDefault( ref.word( "momentumPredictor" ), ref.Switch( True ) )
    
    return pimple, nCorr, nNonOrthCorr, momentumPredictor


#--------------------------------------------------------------------------------------------------------------------------
def setRegionFluidFields( i, fluidRegions, thermoFluid, rhoFluid, kappaFluid, UFluid, \
                          phiFluid, turbulence, KFluid, dpdtFluid, initialMassFluid, ghFluid, ghfFluid, p_rghFluid, radiation ):
    mesh = fluidRegions[ i ]

    thermo = thermoFluid[ i ]
    rho = rhoFluid[ i ]
    kappa = kappaFluid[ i ]
    U = UFluid[ i ]
    phi = phiFluid[ i ]
    
    turb = turbulence[ i ]
    K = KFluid[ i ]
    dpdt = dpdtFluid[ i ]

    p = thermo.p()
    psi = thermo.psi()
    h = thermo.h()
    
    p_rgh = p_rghFluid[ i ]
    gh = ghFluid[ i ]
    ghf = ghfFluid[ i ]
    
    rad = radiation[ i ]
    
    initialMass = ref.dimensionedScalar( ref.word( "massIni" ), ref.dimMass, initialMassFluid[ i ] )

    return mesh, thermo, rho, kappa, K, U, phi, turb, dpdt, p, psi, h, initialMass, p_rgh, gh, ghf, rad


#--------------------------------------------------------------------------------------------------------------------------
def storeOldFluidFields( p, rho ):
    p.storePrevIter()
    rho.storePrevIter()
    pass


#--------------------------------------------------------------------------------------------------------------------------
def initContinuityErrs( fluidRegions ):
    cumulativeContErr = list()
    for index in range( fluidRegions.__len__() ) :
        cumulativeContErr.append( 0.0 )
        pass
    
    return cumulativeContErr


#--------------------------------------------------------------------------------------------------------------------------
def compressibleContinuityErrors( i, mesh, rho, thermo, cumulativeContErr ) :
    totalMass = ref.fvc.domainIntegrate(rho)
       
    sumLocalContErr =  ( ref.fvc.domainIntegrate( ( rho() - thermo.rho()  ).mag() )  / totalMass ).value() # mixed calculations
    
    globalContErr =  ( ref.fvc.domainIntegrate( rho() - thermo.rho() ) / totalMass ).value() # mixed calculations

    cumulativeContErr[i] = cumulativeContErr[i] + globalContErr
    
    ref.ext_Info()<< "time step continuity errors (" << mesh.name() << ")" \
        << ": sum local = " << sumLocalContErr \
        << ", global = " << globalContErr \
        << ", cumulative = " << cumulativeContErr[i] \
        << ref.nl
        
    return cumulativeContErr


#--------------------------------------------------------------------------------------------------------------------------
def fun_UEqn( rho, U, K, phi, ghf, p_rgh, turb, mesh, momentumPredictor, finalIter ):
    # Solve the Momentum equation
    
    UEqn = man.fvm.ddt( rho, U ) + man.fvm.div( phi, U ) + man.fvVectorMatrix( turb.divDevRhoReff( U ), man.Deps( turb, U ) )

    UEqn.relax()
    
    if momentumPredictor :
        ref.solve( UEqn() == ref.fvc.reconstruct( ( -ghf * ref.fvc.snGrad( rho ) - ref.fvc.snGrad( p_rgh ) ) * mesh.magSf() ), 
                   mesh.solver( U.select( finalIter ) ) )
        K << 0.5 * U.magSqr()
        pass
    
    return UEqn


#--------------------------------------------------------------------------------------------------------------------------
def fun_hEqn( rho, h, phi, turb, K, dpdt, thermo, rad,  mesh, oCorr, nOuterCorr, finalIter ):
    
    hEqn = ( ( ref.fvm.ddt( rho, h ) + ref.fvm.div( phi, h ) - ref.fvm.laplacian( turb.alphaEff(), h ) ) \
               == dpdt() - ( ref.fvc.ddt( rho, K ) + ref.fvc.div( phi, K ) ) + rad.Sh( thermo() )  )  # mixed calculation
    
    hEqn.relax()
    hEqn.solve( mesh.solver( h.select( finalIter ) ) ) 
   
    thermo.correct()
    rad.correct()
    
    ref.ext_Info()<< "Min/max T:" << thermo.T().ext_min().value() << ' ' \
        << thermo.T().ext_max().value() << ref.nl
        
    pass


#---------------------------------------------------------------------------------------------------------    
def fun_pEqn( i, mesh, p, rho, turb, thermo, thermoFluid, kappa, K, UEqn, U, phi, psi, dpdt, initialMass, p_rgh, gh, ghf, \
              nNonOrthCorr, oCorr, nOuterCorr, corr, nCorr, cumulativeContErr ) :
    
    closedVolume = p_rgh.needReference()
    
    compressibility = ref.fvc.domainIntegrate( psi )
    
    compressible = ( compressibility.value() > ref.SMALL )

    rho << thermo.rho()
    
    rUA = 1.0 / UEqn.A()
    
    rhorUAf = ref.surfaceScalarField( ref.word( "(rho*(1|A(U)))" ) , ref.fvc.interpolate( rho * rUA ) )

    U << rUA * UEqn.H() 

    phiU = ( ref.fvc.interpolate( rho ) *
                 (  ( ref.fvc.interpolate( U ) & mesh.Sf() ) +
                      ref.fvc.ddtPhiCorr( rUA, rho, U, phi ) ) )
    phi << phiU - rhorUAf * ghf * ref.fvc.snGrad( rho ) * mesh.magSf()

    p_rghDDtEqn = ref.fvc.ddt( rho ) + psi * ref.correction( ref.fvm.ddt( p_rgh ) ) + ref.fvc.div( phi )
    
    # Thermodynamic density needs to be updated by psi*d(p) after the
    # pressure solution - done in 2 parts. Part 1:
    thermo.rho() << thermo.rho() - psi * p_rgh

    for nonOrth in range ( nNonOrthCorr + 1 ):
        p_rghEqn = p_rghDDtEqn - ref.fvm.laplacian( rhorUAf, p_rgh )
        p_rghEqn.solve( mesh.solver( p_rgh.select( ( oCorr == nOuterCorr-1 and corr == nCorr-1 and nonOrth == nNonOrthCorr ) ) ) )
        
        if nonOrth == nNonOrthCorr :
           phi += p_rghEqn.flux()
           pass
        pass
    
    # Second part of thermodynamic density update
    thermo.rho() << thermo.rho() + psi * p_rgh
    
    # Correct velocity field
    U += rUA * ref.fvc.reconstruct( ( phi() - phiU ) / rhorUAf ) # mixed calculations
    U.correctBoundaryConditions()
    K << 0.5 * U.magSqr()
    
    p << p_rgh + rho * gh

    # Update pressure time derivative
    dpdt << ref.fvc.ddt( p )

    # Solve continuity
    ref.rhoEqn( rho, phi )   
    
    # Update continuity errors
    cumulativeContErr = compressibleContinuityErrors( i, mesh, rho, thermo, cumulativeContErr )
    
    # For closed-volume cases adjust the pressure and density levels
    # to obey overall mass continuity
    if closedVolume and compressible:
       p += ( initialMass - ref.fvc.domainIntegrate( thermo.rho() ) ) / compressibility
       rho << thermo.rho()
       p_rgh << p - rho * gh()
       pass
    #Update thermal conductivity
    kappa << thermoFluid[ i ].Cp() * turb.alphaEff()
        
    return cumulativeContErr


#--------------------------------------------------------------------------------------------------------------------------
def solveFluid( i, mesh, thermo, rad, thermoFluid, rho, kappa, K, U, phi, h, turb, dpdt, p, psi, initialMass, p_rgh, gh, ghf,\
                oCorr, nCorr, nOuterCorr, nNonOrthCorr, momentumPredictor,cumulativeContErr, finalIter ) :
    if finalIter:
        mesh.add( ref.word( "finalIteration" ), True )
        pass

    if oCorr == 0 :
        ref.rhoEqn( rho, phi )
        pass
    
    UEqn = fun_UEqn( rho, U, K, phi, ghf, p_rgh, turb, mesh, momentumPredictor, finalIter )
    fun_hEqn( rho, h, phi, turb, K, dpdt, thermo, rad,  mesh, oCorr, nOuterCorr, finalIter )
    
    # --- PISO loop
    for corr in range( nCorr ):
        cumulativeContErr =  fun_pEqn( i, mesh, p, rho, turb, thermo, thermoFluid, kappa, K, UEqn, U, phi, psi, dpdt, initialMass, p_rgh, gh, ghf,
                                       nNonOrthCorr, oCorr, nOuterCorr, corr, nCorr, cumulativeContErr )
        pass
    
    turb.correct()
    rho << thermo.rho()
    
    if finalIter:
        mesh.remove( ref.word( "finalIteration" ) )
        pass
    
    return cumulativeContErr


#-----------------------------------------------------------------------------------------------------------------
