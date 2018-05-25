using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    public abstract class AbstractLazyCollection<TKey, TElem>
    {
        protected string wordForKey { get; private set; }
        protected string wordForElement { get; private set; }
        protected TKey[] keys { get; private set; }
        protected int count { get; private set; }
        private bool[] initialized;
        private TElem[] elements;
        private Dictionary<TKey,int> elementIndex;

        protected TElem GetElement(int i){
            if (i < 0 || i >= count){
                throw new IndexOutOfRangeException($"Cannot get {i}-th {wordForElement} of {count}");
            }
            if (!initialized[i]){
                elements[i] = InitializeElement(i);
                initialized[i] = true;
            }
            return elements[i];            
        }
        protected IEnumerable<TElem> GetElements(){
            return Enumerable.Range(0, count).Select(i => GetElement(i));
        }
        protected bool ContainsKey(TKey key) {
            return count > 0 && elementIndex.ContainsKey(key);
        }
        protected TElem GetElementByKey(TKey key){
            if (!ContainsKey(key)){
                throw new KeyNotFoundExceptionWithKey<TKey>(key, $"This {this.GetType()} does not contain {wordForElement} with {wordForKey} '{key}'");
            }
            return GetElement(elementIndex[key]);
        }

        protected abstract TElem InitializeElement(int i);

        protected AbstractLazyCollection(string wordForKey, string wordForElement){
            this.wordForKey = wordForKey;
            this.wordForElement = wordForElement;
            this.count = 0;
        }

        protected void Initialize(IEnumerable<TKey> keys){
            this.keys = keys.ToArray();
            this.count = this.keys.Length;
            this.initialized = new bool[count];
            this.elements = new TElem[count];
            this.elementIndex = new Dictionary<TKey, int>(keys.Select( (key, i) => new KeyValuePair<TKey,int>(key, i) ));
        }

    }
}