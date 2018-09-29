
using System;
using System.Collections.Generic;
using System.Linq;
using /*SecStrAnnot2.*/Cif.Raw;
using /*SecStrAnnot2.*/Cif.Tables;

namespace /*SecStrAnnot2.*/Cif
{
    public class CifCategory
    {
        private CifParser parser;
        public string[] ItemKeywordNames { get; private set; }
        public CifItem[] Items { get; private set; }
        private Dictionary<string,int> itemIndex;
        public int RowCount { get; private set; }

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
            int[] counts = Items.Select(it => it.Count).Distinct().ToArray();
            if (counts.Length == 1) { // OK
                RowCount = counts[0];
            } else if (counts.Length == 0){ // no items
                RowCount = 0;
            } else { // error
                throw new CifException("Items in category " + Name + " do not have the same number of rows");
            }
        }

        public Table MakeTable(params ValueTuple<string,CifValueType>[] schema) => new Table(this, schema);

        public Table MakeTable(int[] rows, params ValueTuple<string,CifValueType>[] schema) => new Table(this, rows, schema);

        /*public int[] GetRowsGroupedByValues(params string[] itemKeywordNames){
            int[] indices = Enumerable.Range(0, RowCount).ToArray();
            int[][] startsOfGroups = new int[itemKeywordNames.Length][];


            CifItem entity = GetItem("label_entity_id");
            CifItem asym = GetItem("label_asym_id");

            int[] startsOfEntities;
            indices = parser.GroupByValuesInEachRegion(entity.iTag, indices, new int[]{0, entity.Count}, out startsOfEntities);

            int[] startsOfAsyms;
            indices = parser.GroupByValuesInEachRegion(asym.iTag, indices, startsOfEntities, out startsOfAsyms);

            // Lib.WriteLineDebug("entity starts: " + startsOfEntities);
            // Lib.WriteLineDebug("asym starts: " + startsOfAsyms);
            return indices;
        }*/
    }
}