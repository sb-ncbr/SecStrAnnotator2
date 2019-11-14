using System;

namespace Cif
{
    public class CifException : Exception
    {
        public CifException(string message) : base(message){ }
    }
}