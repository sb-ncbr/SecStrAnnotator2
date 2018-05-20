using System;
using System.Collections.Generic;

namespace SecStrAnnot2.Cif
{
    public class TableBunch : Dictionary<string, Table>
    {
        /*private TableBunch(CifBlock block){
            var multidict = new Dictionary<string, List<string>>();
            foreach (string key in block.Keys){
                string[] parts = key.Split('.', 2);
                if (parts.Length < 2){
                    throw new CifException("Tags are expected to contain at least one . character: " + key);
                }
                if (!multidict.ContainsKey(parts[0])){
                    multidict[parts[0]] = new List<string>();
                }
                multidict[parts[0]].Add(parts[1]);
            }
            foreach(var cat_tags in multidict){
                string category = cat_tags.Key;
                List<string> tags = cat_tags.Value;
                this[category] = Table.FromCifBlock(block, category, tags);
            }
            //Console.WriteLine(multidict.EnumerateMultidict());
            //Console.WriteLine(multidict["_atom_site"].Enumerate() );
        }
        public static TableBunch FromCifBlock(CifBlock block){
            return new TableBunch(block);
        }*/
    }
}