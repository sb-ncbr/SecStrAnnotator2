using System;
using System.Collections.Generic;
using System.Linq;

namespace Cif
{
    public class MaybeValueArray<T> where T:struct
    {
        private T[] array;
        private T?[] arrayNullable;
        /* Create infinite MaybeValueArray with all values null. */
        public MaybeValueArray(){
            this.array = null;
            this.arrayNullable = null;
        }
        /* Create MaybeValueArray with all values not-null. */
        public MaybeValueArray(T[] array){
            this.array = array;
            this.arrayNullable = null;
        }
        /* Create MaybeValueArray with some values possibly null. */
        public MaybeValueArray(T?[] arrayNullable){
            this.array = null;
            this.arrayNullable = arrayNullable;
        }
        public bool IsArray => array != null || arrayNullable != null;
        public bool IsNull => array == null && arrayNullable == null;
        public T? this[int i] => array != null ? array[i] : arrayNullable != null ? arrayNullable[i] : null as T?;
    }
}