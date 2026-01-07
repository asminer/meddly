
/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2009, Iowa State University Research Foundation, Inc.

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this library.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
    This program simply displays info about MEDDLY itself.
*/

#include <iostream>
#include <stdlib.h>
#include "../src/meddly.h"

void showDocs(const char* doc)
{
    using namespace std;
    cout << "        ";
    unsigned col = 8;
    for (unsigned i=0; doc[i]; ) {
        // Print next word
        for (; doc[i]; i++) {
            if (' ' == doc[i]) break;
            if ('\n' == doc[i]) break;
            cout << doc[i];
            ++col;
        }
        // Ready for next line?
        if ('\n' == doc[i]) {
            // Paragraph; add a blank line
            cout << "\n";
            ++i;
            col = 99999;
        }
        if (col > 70) {
            cout << "\n        ";
            col = 8;
        } else {
            cout << " ";
            ++col;
        }
        // consume whitespace
        while (' ' == doc[i]) {
            i++;
        }
    }
    cout << "\n";
}

int main()
{
    using namespace std;
    using namespace MEDDLY;

    builtin_init B(nullptr);
    B.setup();

    cout << "===========================================================================\n";
    cout << "Unary operations, for use in apply()\n";
    cout << "===========================================================================\n\n";
    for (unsigned i=0; i<B.numUnary(); i++) {
        const unary_factory *F = B.getUnary(i);
        if (!F) continue;

        cout << "    " << F->getName() << "\n";
        showDocs(F->getDocs());
        cout << "\n";
        cout << "        Implemented in " << F->getFile() << ".\n";
        cout << "\n\n";
    }

    cout << "===========================================================================\n";
    cout << "Binary operations, for use in apply()\n";
    cout << "===========================================================================\n\n";
    for (unsigned i=0; i<B.numBinary(); i++) {
        const binary_factory *F = B.getBinary(i);
        if (!F) continue;

        cout << "    " << F->getName() << "\n";
        showDocs(F->getDocs());
        cout << "\n";
        cout << "        Implemented in " << F->getFile() << ".\n";
        cout << "\n\n";
    }

    return 0;
}
