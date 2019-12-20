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
    public enum JsonType
    {
        Null, 
        String, 
        Double, 
        Int, 
        Bool, 
        Object, 
        List
    }
}