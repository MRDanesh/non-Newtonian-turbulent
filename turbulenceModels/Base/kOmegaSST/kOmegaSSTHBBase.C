/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kOmegaSSTHBBase.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField>
kOmegaSSTHBBase<TurbulenceModel, BasicTurbulenceModel>::kOmegaSSTHBBase::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField>
kOmegaSSTHBBase<TurbulenceModel, BasicTurbulenceModel>::kOmegaSSTHBBase::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField>
kOmegaSSTHBBase<TurbulenceModel, BasicTurbulenceModel>::kOmegaSSTHBBase::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField>
kOmegaSSTHBBase<TurbulenceModel, BasicTurbulenceModel>::kOmegaSSTHBBase::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class TurbulenceModel, class BasicTurbulenceModel>
void kOmegaSSTHBBase<TurbulenceModel, BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& F2
)
{
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F2*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class TurbulenceModel, class BasicTurbulenceModel>
void kOmegaSSTHBBase<TurbulenceModel, BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))), F23());
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField::Internal>
kOmegaSSTHBBase<TurbulenceModel, BasicTurbulenceModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField::Internal>
kOmegaSSTHBBase<TurbulenceModel, BasicTurbulenceModel>::epsilonByk
(
    const volScalarField::Internal& F1,
    const volScalarField::Internal& F2
) const
{
    return betaStar_*omega_();
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<fvScalarMatrix>
kOmegaSSTHBBase<TurbulenceModel, BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<fvScalarMatrix>
kOmegaSSTHBBase<TurbulenceModel, BasicTurbulenceModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<fvScalarMatrix> kOmegaSSTHBBase<TurbulenceModel, BasicTurbulenceModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceModel, class BasicTurbulenceModel>
kOmegaSSTHBBase<TurbulenceModel, BasicTurbulenceModel>::kOmegaSSTHBBase
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    TurbulenceModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )


{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TurbulenceModel, class BasicTurbulenceModel>
bool kOmegaSSTHBBase<TurbulenceModel, BasicTurbulenceModel>::read()
{
    if (TurbulenceModel::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class TurbulenceModel, class BasicTurbulenceModel>
void kOmegaSSTHBBase<TurbulenceModel, BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    
    
    
    

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;
    const volScalarField& nu = this->nu(); // mean viscosity
    fv::options& fvOptions(fv::options::New(this->mesh_));

    BasicTurbulenceModel::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );


    // ********** Read non-Newtonian *************//
    IOdictionary turbulenceProperties
    (
        IOobject
        (
            "turbulenceProperties",
            this->mesh().time().constant(),
            this->mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar K(turbulenceProperties.lookup("K"));
    dimensionedScalar tau0(turbulenceProperties.lookup("tau0"));
    dimensionedScalar n(turbulenceProperties.lookup("n"));
    dimensionedScalar C_beta(turbulenceProperties.lookup("C_beta"));  // For Newtonian case must be zero
    dimensionedScalar C_x(turbulenceProperties.lookup("C_x"));
    dimensionedScalar C_zeta(turbulenceProperties.lookup("C_zeta"));
    dimensionedScalar C_E1(turbulenceProperties.lookup("C_E1"));
    dimensionedScalar C_E2(turbulenceProperties.lookup("C_E2"));
    // **********************//


    tmp<volTensorField> tgradU = fvc::grad(U);

    volScalarField extraNN = C_beta*betaStar_*omega_*k_/nu;  //extra term for non-Newtonian
    volScalarField S2 = 2*magSqr(symm(tgradU()));       // Non-modified mean strain rate (squared)
    volScalarField newS2 = S2 + extraNN;  // Modified mean strain rate (squared)

    volScalarField::Internal GbyNu(dev(twoSymm(tgradU()())) && tgradU()());
    volScalarField::Internal G(this->GName(), nut()*GbyNu);
    tgradU.clear();

    //********************* NOTE **********************//
    // rho is dimensionless, becasue OpenFOAM just take kinematic viscosity. So in the code mu is actually nu
    //****************************************//


    //****************** non-Newtonian calculations **********************//
    volScalarField gammaDot = sqrt(newS2);  // Modified mean strain rate
    dimensionedScalar SMALL_gammaDot("SMALL_gamma", gammaDot.dimensions(), VSMALL);
    gammaDot = max(gammaDot, SMALL_gammaDot); // stabilization
    volScalarField dmu_dgamma = -tau0 / (gammaDot * gammaDot) + K * (n - 1.0) * pow(gammaDot, n - 2.0);

    volScalarField muNN = dmu_dgamma * extraNN/ gammaDot; 
    volScalarField chiNN = -C_x*muNN * S2;

    volScalarField zetaCoeff = dmu_dgamma * S2 / gammaDot;
    volScalarField zetaNN = C_zeta * fvc::laplacian(zetaCoeff, k_);  // check this

    //volScalarField nut_nn = this->nut();

    scalar FE_n = 0.5 * tanh(8.0 * (n.value() - 0.75)) + 0.5;   //Ok
    scalar C_E = C_E1.value() * FE_n + C_E2.value() * (1.0 - FE_n); //Ok
    //********************

    //Info << "nut range: " << gMin(nut) << " / " << gMax(nut) << endl;

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            alpha()*rho()*gamma
           *min
            (
                GbyNu,
                (c1_/a1_)*betaStar_*omega_()
               *max(a1_*omega_(), b1_*F23()*sqrt(S2()))
            )
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + Qsas(S2(), gamma, beta)
          + omegaSource()
          + fvOptions(alpha, rho, omega_)
          + alpha()*(C_E*gamma/nut)*chiNN
          + alpha()*(C_E*gamma/nut)*zetaNN
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        alpha()*rho()*Pk(G)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilonByk(F1, F23), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
      + alpha() * chiNN //ok
      + alpha() * zetaNN //ok
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(S2, F23);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
