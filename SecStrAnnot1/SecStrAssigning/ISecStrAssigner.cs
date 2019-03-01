using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein.SecStrAssigning
{
    public interface ISecStrAssigner {
        SecStrAssignment GetSecStrAssignment(); 
        String GetDescription();
    }	
}

