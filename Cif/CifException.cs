using System;

namespace SecStrAnnot2.Cif
{
    public class CifException : Exception
    {
        public CifException(string message) : base(message){ }
    }
}