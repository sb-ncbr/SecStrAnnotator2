using System;
using System.Collections.Generic;

namespace Cif
{
    public class KeyNotFoundExceptionWithKey<TKey> : KeyNotFoundException
    {
        public TKey Key { get; private set; }

        public KeyNotFoundExceptionWithKey(TKey key) 
            : base() { 
                Key = key;
            }

        public KeyNotFoundExceptionWithKey(TKey key, string message) 
            : base(message) {
                Key = key;
            }

        public KeyNotFoundExceptionWithKey(TKey key, string message, Exception innerException) 
            : base(message, innerException) {
                Key = key;
            }
    }
}