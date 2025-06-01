    // Momentum predictor

    // non-Newtonian stress

	const volScalarField& omega = mesh.lookupObject<volScalarField>("omega");
	const volScalarField& k = mesh.lookupObject<volScalarField>("k");
    const volScalarField& nu = mesh.lookupObject<volScalarField>("nu");

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField extraNN = C_beta*betaStar_*omega*k/nu;  //extra term for non-Newtonian
    volScalarField S2(2*magSqr(symm(tgradU()))+ extraNN);  // Modified mean strain rate (squared)
    volScalarField gammaDot = sqrt(S2);  // Modified mean strain rate
    dimensionedScalar SMALL_gammaDot("SMALL_gamma", gammaDot.dimensions(), VSMALL);
    gammaDot = gammaDot+SMALL_gammaDot;  
    volScalarField dmu_dgamma = -tau0 / (gammaDot * gammaDot) + K * (n - 1.0) * pow(gammaDot, n - 2.0);
    volScalarField muNN = dmu_dgamma * (C_beta* betaStar_ * k * omega / (nu * gammaDot));  // Check omega_ and k_ later
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
