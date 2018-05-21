
using System;
using System.Collections.Generic;
using System.Linq;
using SecStrAnnot2.Cif.Raw;
using SecStrAnnot2.Cif.Tables;

namespace SecStrAnnot2.Cif
{
    public class CifCategory
    {
        private CifParser parser;
        public string[] ItemKeywordNames { get; private set; }
        public CifItem[] Items { get; private set; }
        private Dictionary<string,int> itemIndex;

        public string Name { get; private set; }
        
        public CifItem this[string itemKeywordName] => GetItem(itemKeywordName);
        public CifItem GetItem(string itemKeywordName) {
            if (!itemIndex.ContainsKey(itemKeywordName)){
                throw new KeyNotFoundExceptionWithKey<string>(itemKeywordName, "This " + this.GetType() + " does not contain keyword '" + itemKeywordName + "'");
            }
            return Items[itemIndex[itemKeywordName]];
        }
        public bool ContainsItem(string itemKeywordName) => itemIndex.ContainsKey(itemKeywordName);

        // Not to be instantiated directly. Use CifBlock.GetCategory().
        internal CifCategory(CifParser parser, string categoryName, int[] iTags){
            this.parser = parser;
            this.Name = categoryName;
            int nTags = iTags.Length;
            this.ItemKeywordNames = new string[nTags];
            this.Items = new CifItem[nTags];
            this.itemIndex = new Dictionary<string, int>(nTags);
            for (int i = 0; i < nTags; i++)
            {
                int iTag = iTags[i];
                string keyword = parser.TagNames[iTag].Split('.', 2)[1];
                ItemKeywordNames[i] = keyword;
                Items[i] = new CifItem(parser, keyword, iTag);
                itemIndex[keyword] = i;
            }
            this.ItemKeywordNames = iTags.Select(iTag => parser.TagNames[iTag].Split('.', 2)[1]).ToArray();
        }

        public Table MakeTable(params ValueTuple<string,CifValueType>[] schema) => new Table(this, schema);

        public Table MakeTable(int[] rows, params ValueTuple<string,CifValueType>[] schema) => new Table(this, rows, schema);
    }
}