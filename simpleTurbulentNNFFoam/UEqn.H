    // Momentum predictor

    // non-Newtonian stress

	const volScalarField& omega = mesh.lookupObject<volScalarField>("omega");
	const volScalarField& k = mesh.lookupObject<volScalarField>("k");
    
    

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    // Strain rate magnitude squared
    volScalarField gammaDot = sqrt(S2);  // Effective shear rate
    dimensionedScalar SMALL_gamma("SMALL_gamma", gammaDot.dimensions(), VSMALL);
    gammaDot = sqrt(S2)+SMALL_gamma;  // Effective shear rate
    // Apparent viscosity from real HB (no regularization)
    
    volScalarField muApp = tau0/gammaDot + K * pow(gammaDot, n - 1);
    // d(mu)/d(gammaDot)
    dimensionedScalar SMALL_value("SMALL_gamma", S2.dimensions(), pow(10,-30));
    volScalarField dmu_dgamma = -tau0 / (gammaDot * gammaDot + SMALL_value) + K * (n - 1.0) * pow(gammaDot, n - 2.0);
    // Compute mu^nn
    volScalarField muNN = dmu_dgamma * (C_beta* betaStar_ * k * omega / (muApp * gammaDot));  
    // Optional: Bound it for stability
    muNN = max(muNN, dimensionedScalar("zero", muNN.dimensions(), 0.0));
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
