using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif
{
    public class Table
    {
        public string Category { get; private set; }
        public int Count{ get; private set; }
        public string[] Keywords{ get; private set; }
        private Dictionary<string, Type> types;
        private Dictionary<string, object> columns;


        public Type GetType(string keyword){
            return types[keyword];
        }
        public T[] GetColumn<T>(string keyword){
            return columns[keyword] as T[];
        }

        private Table (IEnumerable<Tuple<string,Type>> schema, int count){
            Keywords = schema.Select(kt => kt.Item1).ToArray();
            types = new Dictionary<string, Type>();
            columns = new Dictionary<string, object>();
            foreach(var keyword_type in schema){
                string keyword = keyword_type.Item1;
                Type type = keyword_type.Item2;
                types[keyword] = type;
                columns[keyword] = Activator.CreateInstance(type.MakeArrayType(), count);
            }
        }
    }
}