using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Threading.Tasks;

using Cif.Components;
using protein.SSEs;
using protein.Libraries;

namespace protein.Annotating
{
    interface IAnnotator
    {
        List<(int, int)> GetMatching();
    }
}