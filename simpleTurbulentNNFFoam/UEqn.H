    // Momentum predictor

    // non-Newtonian stress

	const volScalarField& omega = mesh.lookupObject<volScalarField>("omega");
	const volScalarField& k = mesh.lookupObject<volScalarField>("k");
    
    Info<< "Readinggggg   1" << endl;

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    Info<< "Readinggggg   2" << endl;
    // Strain rate magnitude squared
    volScalarField gammaDot = sqrt(S2);  // Effective shear rate
    dimensionedScalar SMALL_gamma("SMALL_gamma", gammaDot.dimensions(), VSMALL);
    gammaDot = sqrt(S2)+SMALL_gamma;  // Effective shear rate
    // Apparent viscosity from real HB (no regularization)
    Info<< "Readinggggg   3" << endl;
    dimensionedScalar SMALL_value1("SMALL_gamma1", gammaDot.dimensions(), pow(10,-30));
    volScalarField muApp = tau0/gammaDot + K * pow(gammaDot+SMALL_value1, n - 1);
    // d(mu)/d(gammaDot)
    dimensionedScalar SMALL_value2("SMALL_gamma2", S2.dimensions(), pow(10,-30));
    volScalarField dmu_dgamma = -tau0 / (gammaDot * gammaDot + SMALL_value2) + K * (n - 1.0) * pow(gammaDot+SMALL_value1, n - 2.0);
    // Compute mu^nn
    Info<< "Readinggggg   4" << endl;
    volScalarField muNN = dmu_dgamma * (C_beta* betaStar_ * k * omega / (muApp * gammaDot));  
    // Optional: Bound it for stability
    // muNN = max(muNN, dimensionedScalar("zero", muNN.dimensions(), 0.0));
    tgradU.clear();

    

   // Replace stress term
   volScalarField nuEffCustom = C_tau * muNN;
   //volVectorField nonNewtonianViscousTerm = 2*fvm::laplacian(nu_nn, U);
    //

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + g
      + turbulence->divDevReff(U)
      + 2*fvm::laplacian(nuEffCustom, U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }
