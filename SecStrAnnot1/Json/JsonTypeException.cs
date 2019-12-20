using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Runtime.Serialization;
using System.Text;
using System.IO;

namespace protein.Json
{
    class JsonTypeException : Exception
    {
        public JsonTypeException(String msg) : base(msg) { }
    }
}