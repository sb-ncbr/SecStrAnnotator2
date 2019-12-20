using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

using Cif.Components;
using Cif.Tables;
using protein.Geometry;

namespace protein.HydrogenAdding
{
    public interface IHydrogenAdder
    {
        // List<Residue> AddHydrogens (IEnumerable<Residue> residues);
        Protein AddHydrogens(Protein protein);
    }
}

