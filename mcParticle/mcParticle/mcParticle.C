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

\*---------------------------------------------------------------------------*/

#include "mcParticle.H"

#include "mcParticleCloud.H"
#include "OStringStream.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(mcParticle, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
//- 2020.09.05@Zmeng
Foam::mcParticle::trackingData::trackingData(mcParticleCloud& mcpc, scalar trackTime)
//Foam::mcParticle::trackData::trackData(mcParticleCloud& mcpc, scalar trackTime)
:
    base(mcpc),
    cloud_(mcpc)
{}


Foam::mcParticle::mcParticle
(
    const mcParticleCloud& c,
    const vector& position,
    const label   celli,
    const scalar  m,
    const vector& UParticle,
    const scalarField&  Phi,
    const vector& shift,
    const label   ghost
)
:
    base(c.pMesh(), position, celli),
    m_(m),
    UParticle_(UParticle),
    Ucorrection_(vector::zero),
    Utracking_(UParticle),
    Omega_(0.0),
    rho_(0.0),
    eta_(1.0),
    shift_(shift),
    Co_(0.0),
    reflectionBoundaryVelocity_(vector::zero),
    ghost_(ghost),
    nSteps_(0),
    isOnInletBoundary_(false),
    reflectedAtOpenBoundary_(false),
    Phi_(Phi)
{
    const polyMesh& mesh = c.mesh();
    meshTools::constrainDirection(mesh, mesh.geometricD(), Utracking_);
    c.computeCourantNo(*this);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//- 2020.09.05@Zmeng
bool Foam::mcParticle::move(mcParticleCloud& cloud, mcParticle::trackingData& td, const scalar trackTime)
//bool Foam::mcParticle::move(mcParticle::trackData& td, const scalar trackTime)
{

    td.switchProcessor = false;
    td.keepParticle = true;

    mcParticleCloud& mcpc = td.cloud();

    if (isOnInletBoundary_)
    {
        stepFraction() = mcpc.random().scalar01();
    }

    const polyMesh& mesh = mcpc.pMesh();
//- 2020.09.05@Zmeng
//    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    scalar tEnd = (1.0 - stepFraction())*trackTime;
//- 2020.09.05@Zmeng
//    scalar dtMax = tEnd;

    isOnInletBoundary_ = false;

//- 2020.09.05@Zmeng IMPORTANT

    while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
    {
        if (debug)
        {
            Pout<< "Time = " << mesh.time().timeName()
                << "  trackTime = " << trackTime
                << "  tEnd = " << tEnd
                << "  stepFraction() = " << stepFraction()
                << "  origId() = " << origId()
                << "  position() = " << position();
        }

        const scalar f = 1 - stepFraction();
        trackToAndHitFace(f*trackTime*Utracking_, f, cloud, td);

        ++nSteps_;

        // if we made too many very small steps, drop the particle
        if (nSteps_ > 1000)
        {
#ifdef FULLDEBUG
            Perr<< "DEBUG: particle " << origId_ << " made more than 1000 "
                   "steps, droping it. Info:\n"
                << info() << nl;
#endif
            mcpc.notifyLostParticle(*this);
            td.keepParticle = false;
            break;
        };
    }

//////////////////////////////////////////////////
/*
    while (td.keepParticle && !td.switchProcessor && tEnd > 0)
    {
        if (debug)
        {
            Pout<< "Time = " << mesh.time().timeName()
                << "  trackTime = " << trackTime
                << "  tEnd = " << tEnd
                << "  stepFraction() = " << stepFraction()
                << "  origId() = " << origId()
                << "  position() = " << position();
        }

        // set the lagrangian time-step
        scalar dt = min(dtMax, tEnd);
        point destPos = position() + dt * Utracking_;
        meshTools::constrainDirection(mesh, mesh.geometricD(), destPos);
        scalar tf = trackToFace(destPos, td);
        ++nSteps_;

        // if we made too many very small steps, drop the particle
        if (nSteps_ > 1000)
        {
#ifdef FULLDEBUG
            Perr<< "DEBUG: particle " << origId_ << " made more than 1000 "
                   "steps, droping it. Info:\n"
                << info() << nl;
#endif
            mcpc.notifyLostParticle(*this);
            td.keepParticle = false;
            break;
        };

        dt *= tf;
        if (debug)
        {
            Pout<< "  trackFraction = " << tf
                << endl;
        }

        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/trackTime;

        if (onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[patch(face())]))
            {
                td.switchProcessor = true;
            }
        }
    }
*/
    return td.keepParticle;
}

//- 2020.09.05@Zmeng
bool Foam::mcParticle::hitPatch(mcParticleCloud&, trackingData&)
{
    return false;
}

// Pre-action before hitting patches
bool Foam::mcParticle::hitPatch
(
    const polyPatch&           patch,
//- 2020.09.05@Zmeng
    mcParticle::trackingData&     td,
//    mcParticle::trackData&     td,
    const label                patchI,
    const scalar               trackFraction,
    const tetIndices&          tetIs
)
{
    if (isA<wedgePolyPatch>(patch))
    {
        const polyMesh& mesh = td.cloud().pMesh();
//- 2020.09.05@Zmeng IMPORTANT
        vector tmpPosition = position();
        meshTools::constrainDirection(mesh, mesh.geometricD(), tmpPosition);
//        meshTools::constrainDirection(mesh, mesh.geometricD(), position());
    }
    else
    {
        mcParticleCloud& c = const_cast<mcParticleCloud&>
        (
            refCast<const mcParticleCloud>(td.cloud())
        );
        c.hitPatch(*this, td, patchI, trackFraction, tetIs);
    }
    return true;
}


void Foam::mcParticle::transformProperties (const tensor& T)
{
    base::transformProperties(T);
    // Only transform fluctuating velocity
    UParticle_ = transform(T, UParticle_);
    Ucorrection_ = transform(T, Ucorrection_);
    Utracking_ = transform(T, Utracking_);
}


void Foam::mcParticle::transformProperties(const vector& separation)
{
    base::transformProperties(separation);
}


Foam::string Foam::mcParticle::info() const
{
    OStringStream oss;
    oss << "Particle Id: " << origId() << ": "
        << "X     = " << position() << ", "
        << "cell  = " << cell() << ", "
        << "m     = " << m() << nl
        << "Ucorrection = " << Ucorrection() << ", "
        << "Utracking = " << Utracking() << ", "
        << "U     = " << UParticle()  << ", "
        << "Phi   = " << Phi() << ", "
        << "ghost = " << ghost() << ", "
        << "shift = " << shift()
        << endl;
    return oss.str();
}


// ************************************************************************* //
