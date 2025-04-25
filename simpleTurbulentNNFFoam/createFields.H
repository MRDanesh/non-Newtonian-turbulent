Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// ********** Read non-Newtonian *************//


    IOdictionary turbulenceProperties
(
    IOobject
    (
        "turbulenceProperties",						
        runTime.constant(),							
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

    dimensionedScalar K(turbulenceProperties.lookup("K"));
    dimensionedScalar tau0(turbulenceProperties.lookup("tau0"));
    dimensionedScalar n(turbulenceProperties.lookup("n"));
    dimensionedScalar betaStar_(turbulenceProperties.lookup("betaStar_"));
    dimensionedScalar C_beta(turbulenceProperties.lookup("C_beta"));
    dimensionedScalar C_tau(turbulenceProperties.lookup("C_tau"));

// Gravity 

#include "readGravitationalAcceleration.H"

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);



#include "createMRF.H"
#include "createFvOptions.H"
