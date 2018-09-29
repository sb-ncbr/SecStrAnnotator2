
using System;
using System.Collections.Generic;
using System.Linq;
using /*SecStrAnnot2.*/Cif.Raw;

namespace /*SecStrAnnot2.*/Cif
{
    public class CifBlock
    {
        private CifParser parser;
        private int iBlock;
        public string[] CategoryNames { get; private set; }
        public string[] ItemFullNames { get; private set; }
        public CifCategory[] Categories { get; private set; }
        public CifItem[] Items { get; private set; }
        private Dictionary<string,int> categoryIndex;
        private Dictionary<string,int> itemIndex;
        private Dictionary<string,string[]> keywordsByCategory;


        public string Name => parser.BlockNames[iBlock];

        public CifCategory this[string categoryName] => GetCategory(categoryName);
        public CifCategory GetCategory(string categoryName) {
            if (!categoryIndex.ContainsKey(categoryName)){
                throw new KeyNotFoundExceptionWithKey<string>(categoryName, "This " + this.GetType() + " does not contain category with name '" + categoryName + "'");
            }
            return Categories[categoryIndex[categoryName]];
        }
        public bool ContainsCategory(string categoryName) => categoryIndex.ContainsKey(categoryName);

        public CifItem GetItem(string itemFullName) {
            if (!itemIndex.ContainsKey(itemFullName)){
                throw new KeyNotFoundExceptionWithKey<string>(itemFullName, "This " + this.GetType() + " does not contain item with full name '" + itemFullName + "'");
            }
            return Items[itemIndex[itemFullName]];
        }
        public bool ContainsItem(string itemFullName) => itemIndex.ContainsKey(itemFullName);

        // Not to be instantiated directly. Use CifPackage.GetBlock().
        internal CifBlock(CifParser parser, int iBlock){
            this.parser = parser;
            this.iBlock = iBlock;
            this.ItemFullNames = parser.TagNames;
            
            // sort items by category
            string[] categoryNames;
            Dictionary<string, string[]> keywordsByCategory_;
            Dictionary<string, int[]> tagIndicesByCategory;
            SortTagsToCategories(ItemFullNames, out categoryNames, out keywordsByCategory_, out tagIndicesByCategory);
            this.CategoryNames = categoryNames;
            this.keywordsByCategory = keywordsByCategory_;

            // initialize categories
            this.Categories = new CifCategory[CategoryNames.Length];
            this.categoryIndex = new Dictionary<string, int>(CategoryNames.Length);
            for (int i = 0; i < CategoryNames.Length; i++)
            {
                string name = CategoryNames[i];
                categoryIndex[name] = i;
                Categories[i] = new CifCategory(parser, name, tagIndicesByCategory[name]);
            }

            // initialize items
            this.Items = new CifItem[ItemFullNames.Length];
            this.itemIndex = new Dictionary<string,int>(ItemFullNames.Length);
            for (int j = 0; j < ItemFullNames.Length; j++)
            {
                string name = ItemFullNames[j];
                itemIndex[name] = j;
                string[] parts = name.Split('.', 2);
                Items[j] = GetCategory(parts[0]).GetItem(parts[1]);
            }
        }

        private static void SortTagsToCategories(string[] tags, 
            out string[] categories, 
            out Dictionary<string,string[]> keywordsByCategory, 
            out Dictionary<string,int[]> tagIndicesByCategory)
        {
            var categoryList = new List<string>();
            var multidictKeywords = new Dictionary<string, List<string>>();
            var multidictIndices = new Dictionary<string, List<int>>();
            
            for (int i = 0; i < tags.Length; i++)
            {
                string tag = tags[i];
                string[] parts = tag.Split('.', 2);
                if (parts.Length < 2){
                    throw new CifException("Tags are expected to have form category.keyword: " + tag);
                }
                string category = parts[0];
                string keyword = parts[1];
                if (!multidictKeywords.ContainsKey(category)){
                    categoryList.Add(category);
                    multidictKeywords[category] = new List<string>();
                    multidictIndices[category] = new List<int>();
                }
                multidictKeywords[category].Add(keyword);
                multidictIndices[category].Add(i);
            }
            categories = categoryList.ToArray();
            keywordsByCategory = new Dictionary<string,string[]>(multidictKeywords.Select(kv => new KeyValuePair<string,string[]>(kv.Key, kv.Value.ToArray())));
            tagIndicesByCategory = new Dictionary<string,int[]>(multidictIndices.Select(kv => new KeyValuePair<string,int[]>(kv.Key, kv.Value.ToArray())));
        }
    }
}