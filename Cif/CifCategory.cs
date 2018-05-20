
using System;
using System.Collections.Generic;
using System.Linq;
using SecStrAnnot2.Cif.Raw;

namespace SecStrAnnot2.Cif
{
    public class CifCategory
    {
        private CifParser parser;
        public string[] ItemKeywordNames { get; private set; }
        public MyCifItem[] Items { get; private set; }
        private Dictionary<string,int> itemIndex;

        public string Name { get; private set; }
        
        public MyCifItem this[string itemKeywordName] => GetItem(itemKeywordName);
        public MyCifItem GetItem(string itemKeywordName) {
            if (!itemIndex.ContainsKey(itemKeywordName)){
                throw new KeyNotFoundException("This " + this.GetType() + " does not contain keyword '" + itemKeywordName + "'");
            }
            return Items[itemIndex[itemKeywordName]];
        }

        internal CifCategory(CifParser parser, string categoryName, int[] iTags){
            this.parser = parser;
            this.Name = categoryName;
            int nTags = iTags.Length;
            this.ItemKeywordNames = new string[nTags];
            this.Items = new MyCifItem[nTags];
            this.itemIndex = new Dictionary<string, int>(nTags);
            for (int i = 0; i < nTags; i++)
            {
                int iTag = iTags[i];
                string keyword = parser.TagNames[iTag].Split('.', 2)[1];
                ItemKeywordNames[i] = keyword;
                Items[i] = new MyCifItem(parser, keyword, iTag);
                itemIndex[keyword] = i;
            }
            this.ItemKeywordNames = iTags.Select(iTag => parser.TagNames[iTag].Split('.', 2)[1]).ToArray();
        }
    }
}