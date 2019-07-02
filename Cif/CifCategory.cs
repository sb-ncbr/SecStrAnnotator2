
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Cif.Raw;
using Cif.Tables;

namespace Cif
{
    public class CifCategory
    {
        internal CifParser parser;
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

        public Table MakeTable(params (string, CifValueType)[] schema) => new Table(this, schema);

        public Table MakeTable(int[] rows, params (string, CifValueType)[] schema) => new Table(this, rows, schema);

        public string MakeCifString() => MakeCifString(Enumerable.Range(0, RowCount).ToArray());
        public string MakeCifString(int[] rows){
            string CATEGORY_SEPARATOR = "#\n";
            string COMMENT_STRING = "# ";
            string LOOP_STRING = "loop_\n";
            string SEPARATOR = " ";
            string ROW_SEPARATOR = "\n";

            StringBuilder builder = new StringBuilder();
            builder.Append(CATEGORY_SEPARATOR);

            if (rows.Length == 0) {
                foreach(CifItem item in this.Items){
                    builder.Append(COMMENT_STRING);
                    builder.Append(item.FullName);
                    builder.Append(ROW_SEPARATOR);
                }
            } else if (rows.Length == 1) {
                // print in 1-value format (without loop)
                foreach(CifItem item in this.Items){
                    builder.Append(item.FullName);
                    builder.Append(SEPARATOR);
                    builder.Append(item.GetStrings()[0]);
                    builder.Append(ROW_SEPARATOR);
                }
            } else {
                // print in multi-value format (loop)
                builder.Append(LOOP_STRING);
                foreach(CifItem item in Items){
                    builder.Append(item.FullName);
                    builder.Append(ROW_SEPARATOR);
                }
                //Version with printing separate values:
                // foreach(int iRow in rows){
                //     foreach(CifItem item in Items){
                //         (int start, int stop) = parser.GetValuePosition(item.iTag, iRow);
                //         builder.Append(parser.Text, start, stop - start);
                //         builder.Append(SEPARATOR);
                //     }
                //     builder.Append(ROW_SEPARATOR);
                // }
                //Version with copying rows directly from the original text:
                int firstTag = Items[0].iTag;
                int lastTag = Items[Items.Length-1].iTag;
                foreach(int iRow in rows){
                    int start = parser.GetValuePosition(firstTag, iRow).Item1;
                    int stop = parser.GetValuePosition(lastTag, iRow).Item2;
                    builder.Append(parser.Text, start, stop - start);
                    builder.Append(ROW_SEPARATOR);
                }
            }
            builder.Append(ROW_SEPARATOR);
            return builder.ToString();
        }

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