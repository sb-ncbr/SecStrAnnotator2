using System;
using System.Collections.Generic;
using System.Linq;

namespace Cif
{
    public class MaybeClassArray<T> where T:class
    {
        private T[] array;
        public MaybeClassArray(T[] array){
            this.array = array;
        }
        public bool IsArray => array != null;
        public bool IsNull => array == null;
        public T this[int i] => IsArray ? array[i] : null;
    }
}