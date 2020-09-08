Info<< "Reading thermophysical properties\n" << endl;

//- 2020.09.07@Zmeng

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);
fluidThermo& fThermo = pThermo();
fThermo.validate(args.executable(), "h", "e");


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

    Info<< "Reading thermophysical properties\n" << endl;
    mcThermo thermo(mesh);
    thermo.correct();

    volScalarField& rho = thermo.rho();
    volScalarField& p = thermo.p();

    #include "compressibleCreatePhi.H"

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

    Info<< "Creating turbulence model\n" << endl;


/*
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );
*/
Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        pThermo
    )
);

    thermo.createCloud();

    dimensionedScalar initialMass = fvc::domainIntegrate(rho);





















