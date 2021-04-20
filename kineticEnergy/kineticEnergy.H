/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 
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

Class
    Foam::functionObjects::kineticEnergy

Group

Description
    Calculates the total kinetic energy of the internal field and writes to a file.
    \verbatim
    kineticEnergy1
    {
        type           kineticEnergy;
        libs ("libkineticEnergyFunctionObject.so");
        ...
        wordData       someWord;
        scalarData     1.0;
        labelData      1;
    }
    \endverbatim

Usage
    \table
        Property     | Description               | Required | Default value
        type         | type name: kineticEnergy | yes      |
    \endtable

SourceFiles
    kineticEnergy.C

\*---------------------------------------------------------------------------*/

#ifndef kineticEnergy_H
#define kineticEnergy_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class kineticEnergy Declaration
\*---------------------------------------------------------------------------*/

class kineticEnergy
:
    public fvMeshFunctionObject,
    public writeFile
{
public:
        virtual bool calc();
protected:

    // Protected Member Functions

        //- File header information
    virtual void writeFileHeader(Ostream& os) const;
    virtual scalar calc(Ostream& os) const;




public:


    //- Runtime type information
    TypeName("kineticEnergy");


    // Constructors

        //- Construct from Time and dictionary
        kineticEnergy
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~kineticEnergy();


    // Member Functions

        //- Read the kineticEnergy data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the kineticEnergy
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
