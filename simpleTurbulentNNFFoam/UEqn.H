    // Momentum predictor

    // non-Newtonian stress

	const volScalarField& omega = mesh.lookupObject<volScalarField>("omega");
	const volScalarField& k = mesh.lookupObject<volScalarField>("k");
    const volScalarField& nu = mesh.lookupObject<volScalarField>("nu");

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField extraNN = C_beta*betaStar_*omega*k/nu;  //extra term for non-Newtonian
    volScalarField S2 = 2*magSqr(symm(tgradU()));       // Non-modified mean strain rate (squared)
    volScalarField newS2 = S2 + extraNN;  // Modified mean strain rate (squared)
    volScalarField gammaDot = sqrt(newS2);  // Modified mean strain rate
    dimensionedScalar SMALL_gammaDot("SMALL_gamma", gammaDot.dimensions(), VSMALL);
    gammaDot = max(gammaDot, SMALL_gammaDot); // stabilization
    volScalarField dmu_dgamma = -tau0 / (gammaDot * gammaDot) + K * (n - 1.0) * pow(gammaDot, n - 2.0);

    volScalarField muNN = dmu_dgamma * (C_beta*betaStar_*omega*k)/ (nu*gammaDot); 
        

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + g
      + turbulence->divDevReff(U)
      + fvc::div(2.0*C_tau*muNN *symm(tgradU()))
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    tgradU.clear();
    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }
