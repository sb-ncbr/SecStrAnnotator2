using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Runtime.Serialization;
using System.Text;
using System.IO;

namespace protein.Json.Helpers
{
    class JsonEnumerator : IEnumerator<JsonValue>
    {
        private IEnumerator<object> inner;

        public JsonEnumerator(IEnumerator<object> innerEnumerator)
        {
            inner = innerEnumerator;
        }

        public JsonValue Current { get { return new JsonValue(inner.Current); } }
        object System.Collections.IEnumerator.Current { get { return Current; } }
        public void Dispose() { inner.Dispose(); }
        public bool MoveNext() { return inner.MoveNext(); }
        public void Reset() { inner.Reset(); }
    }

}