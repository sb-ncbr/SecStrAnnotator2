using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    public class Structure
    {
        private const string CATEGORY_NAME = "_atom_site";

        public AtomTable Atoms { get; private set; }
        public ResidueTable Residues { get; private set; }
        public FragmentTable Fragments { get; private set; }
        public ChainTable Chains { get; private set; }
        public EntityTable Entities { get; private set; }

        public static Structure FromCifBlock(CifBlock cifBlock){
            if (!cifBlock.ContainsCategory(CATEGORY_NAME)){
                throw new CifException("Cannot read a structure from a CIF block " + cifBlock.Name + ", because it does not contain category " + CATEGORY_NAME);
            }
            return FromCifCategory(cifBlock[CATEGORY_NAME]);
        }
        public static Structure FromCifBlock(CifBlock cifBlock, int[] rows){
            if (!cifBlock.ContainsCategory(CATEGORY_NAME)){
                throw new CifException("Cannot read a structure from a CIF block " + cifBlock.Name + ", because it does not contain category " + CATEGORY_NAME);
            }
            return FromCifCategory(cifBlock[CATEGORY_NAME], rows);
        }

        public static Structure FromCifCategory(CifCategory cifCategory){
            int[] rows = Enumerable.Range(0, cifCategory.RowCount).ToArray();
            return FromCifCategory(cifCategory, rows);
        }
        public static Structure FromCifCategory(CifCategory cifCategory, int[] rows){
            return new Structure(cifCategory, rows);
        }

        private Structure(CifCategory cifCategory, int[] rows){
            throw new NotImplementedException();
        }
    }
}