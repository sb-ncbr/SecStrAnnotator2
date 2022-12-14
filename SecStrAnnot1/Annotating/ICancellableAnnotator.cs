using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Threading.Tasks;

using Cif.Components;
using protein.Sses;
using protein.Libraries;

namespace protein.Annotating
{
    interface ICancellableAnnotator: IAnnotator
    {
        List<(int, int)> GetMatching(System.Threading.CancellationToken cancellationToken);
    }
}