using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein.SecStrAssigning.Helpers
{
    interface IHBondFinder{
        bool IsHBond (int donor, int acceptor);
        List<int> FindHAcceptors(int donor);
        List<int> FindHDonors(int acceptor);
    }	
}

