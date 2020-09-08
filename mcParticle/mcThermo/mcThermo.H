/*---------------------------------------------------------------------------*\
                pdfFoam: General Purpose PDF Solution Algorithm
                   for Reactive Flow Simulations in OpenFOAM

 Copyright (C) 2012 Michael Wild, Heng Xiao, Patrick Jenny,
                    Institute of Fluid Dynamics, ETH Zurich
-------------------------------------------------------------------------------
License
    This file is part of pdfFoam.

    pdfFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) version 3 of the same License.

    pdfFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with pdfFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::mcThermo

Description
    Thermo class based on MC particles

Note
    The constructor does @b not create the particle cloud, because
    mcParticleCloud requires a turbulence model, but the compressible
    turbulence models also required mcThermo to be passed in their constructor.
    To break this dependency loop, it is the users responsibility to call
    mcThermo::createCloud() after creating the turbulence model.

SourceFiles
    mcThermo.C

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#ifndef mcThermo_H
#define mcThermo_H

#include "basicThermo.H"
#include "mcParticleCloud.H"
#include "volFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class mcThermo Declaration
\*---------------------------------------------------------------------------*/

class mcThermo
:
    public basicThermo
{
    // Private data

        //- The mesh
        const fvMesh& mesh_;

        //- Particle cloud
        autoPtr<mcParticleCloud> cloudP_;

        //- Density field
        volScalarField rho_;

        //- Viscosity
        dimensionedScalar nu_;

    // Private member functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        mcThermo(const mcThermo&);


public:

    //- Runtime type information
    TypeName("mcThermo");


    // Constructors

        //- Construct from mesh
        mcThermo(const fvMesh& mesh, const word& phaseName=word::null);


    //- Destructor
    virtual ~mcThermo();


    // Member functions

        //- Create the particle cloud
        void createCloud();

        //- Evolve the cloud
        // @returns The maximum residual
        scalar evolve();

        //- Update properties
        virtual void correct();

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();



//- 2020.09.05@Zmeng

            //- Return the name of the thermo physics
            virtual word thermoName() const
            {
                NotImplemented;
                return word::null;
            }
            //- Return true if the equation of state is incompressible
            //  i.e. rho != f(p)
            virtual bool incompressible() const
            {
                NotImplemented;
                return true;
            }

            //- Return true if the equation of state is isochoric
            //  i.e. rho = const
            virtual bool isochoric() const
            {
                NotImplemented;
                return false;
            }

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he()
            {
                NotImplemented;
                
            }

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const
            {
                NotImplemented;
                return volScalarField::null();
            }

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const
            {
                NotImplemented;
                return tmp<volScalarField>();
            }

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const
            {
                NotImplemented;
                return tmp<volScalarField>();
            }

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const
            {
                NotImplemented;
                return tmp<volScalarField>();
            }

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const
            {
                NotImplemented;
                return tmp<volScalarField>();
            }

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }

            //- Gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const
            {
                NotImplemented;
                return tmp<volScalarField>();
            }

            //- Gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> Cpv() const
            {
                NotImplemented;
                return tmp<volScalarField>();
            }

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }

            //- Heat capacity ratio []
            virtual tmp<volScalarField> CpByCpv() const
            {
                NotImplemented;
                return tmp<volScalarField>();
            }

            //- Heat capacity ratio for patch []
            virtual tmp<scalarField> CpByCpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const
            {
                NotImplemented;
                return tmp<volScalarField>();
            }

            //- Molecular weight for patch [kg/kmol]
            virtual tmp<scalarField> W(const label patchi) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }


        // Access to transport state variables

            //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
            virtual const volScalarField& alpha() const
            {
                NotImplemented;
                return volScalarField::null();
            }

            //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
            virtual const scalarField& alpha
            (
                const label patchi
            ) const
            {
                NotImplemented;
                return volScalarField::null();
            }


        // Fields derived from transport state variables

            //- Thermal diffusivity for temperature of mixture [W/m/K]
            virtual tmp<volScalarField> kappa() const
            {
                NotImplemented;
                return tmp<volScalarField>();
            }

            //- Thermal diffusivity for temperature of mixture
            //  for patch [W/m/K]
            virtual tmp<scalarField> kappa
            (
                const label patchi
            ) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }

            //- Thermal diffusivity for energy of mixture [kg/m/s]
            virtual tmp<volScalarField> alphahe() const
            {
                NotImplemented;
                return tmp<volScalarField>();
            }

            //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphahe(const label patchi) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture [W/m/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField&
            ) const
            {
                NotImplemented;
                return tmp<volScalarField>();
            }

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }

            //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const
            {
                NotImplemented;
                return tmp<volScalarField>();
            }

            //- Effective thermal turbulent diffusivity of mixture
            //  for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const
            {
                NotImplemented;
                return tmp<scalarField>();
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
