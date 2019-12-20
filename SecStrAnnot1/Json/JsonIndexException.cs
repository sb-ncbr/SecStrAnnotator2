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
    class JsonIndexException : Exception
    {
        public int Index { get; private set; }
        public int Count { get; private set; }

        public JsonIndexException(String msg, int index, int count) : base(msg)
        {
            Index = index;
            Count = count;
        }
    }
}