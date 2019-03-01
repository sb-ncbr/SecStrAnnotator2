using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein.SecStrAssigning
{
    public class SecStrAssignmentException : Exception{
        public SecStrAssignmentException(String message) : base(message) {}
        public SecStrAssignmentException(String message, Exception innerException) : base(message,innerException) {}
    }	
}

