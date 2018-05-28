using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    public class ModelCollection : AbstractLazyCollection<int, Model>
    {
        public const string CATEGORY_NAME = "_atom_site";

        public int[] ModelNumbers => base.keys;
        public int Count => base.count;
        public Model GetModel(int i) => base.GetElement(i);
        public IEnumerable<Model> GetModels() => base.GetElements();
        public bool ContainsModel(int modelNumber) => base.ContainsKey(modelNumber);
        public Model GetModelByModelNumber(int modelNumber) => base.GetElementByKey(modelNumber);

        private CifCategory category;
        private int[] rows;
        private int[] modelStartRowIndex;

        protected override Model InitializeElement(int i){
            int rowCount = modelStartRowIndex[i+1] - modelStartRowIndex[i];
            int[] modelRows = new int[rowCount];
            Array.Copy(rows, modelStartRowIndex[i], modelRows, 0, rowCount);
            return new Model(category, modelRows, ModelNumbers[i]);
        }

        private ModelCollection(CifCategory category, int[] rows) : base("model number", "model") {
            this.category = category;
            this.rows = rows;
            if (category.ContainsItem(Model.MODEL_NUM_COLUMN)) {
                CifItem modelNumItem = category[Model.MODEL_NUM_COLUMN];
                int[] modelNumberPerAtom = modelNumItem.GetIntegers(rows);
                List<int> runStarts;
                try {
                    runStarts = Lib.RunStartsInOrderedArray(modelNumberPerAtom);
                } catch (ArgumentException e) {
                    throw new NotImplementedException("Model numbers (" + modelNumItem.FullName + ") are not in increasing order (not supported by current version)", e);
                }
                IEnumerable<int> modelNumbers = runStarts.Select(start => modelNumberPerAtom[start]);
                modelStartRowIndex = runStarts.Append(rows.Length).ToArray();
                base.Initialize(modelNumbers);
            } else if (rows.Length > 0) {
                modelStartRowIndex = new int[]{ 0, rows.Length };
                base.Initialize(new int[] { Model.DEFAULT_MODEL_NUM });
            } else {
                modelStartRowIndex = new int[] { 0 };
                base.Initialize(new int[]{});
            }
        }

        public static ModelCollection FromCifBlock(CifBlock cifBlock){
            if (!cifBlock.ContainsCategory(CATEGORY_NAME)){
                throw new CifException("Cannot read a structure from a CIF block " + cifBlock.Name + ", because it does not contain category " + CATEGORY_NAME);
            }
            return FromCifCategory(cifBlock[CATEGORY_NAME]);
        }
        public static ModelCollection FromCifBlock(CifBlock cifBlock, int[] rows){
            if (!cifBlock.ContainsCategory(CATEGORY_NAME)){
                throw new CifException("Cannot read a structure from a CIF block " + cifBlock.Name + ", because it does not contain category " + CATEGORY_NAME);
            }
            return FromCifCategory(cifBlock[CATEGORY_NAME], rows);
        }

        public static ModelCollection FromCifCategory(CifCategory cifCategory){
            int[] rows = Enumerable.Range(0, cifCategory.RowCount).ToArray();
            return FromCifCategory(cifCategory, rows);
        }
        public static ModelCollection FromCifCategory(CifCategory cifCategory, int[] rows){
            return new ModelCollection(cifCategory, rows);
        }
    }
}