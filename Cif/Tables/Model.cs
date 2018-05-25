using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    public class Model
    {
        public int ModelNumber { get; private set; }
        public AtomTable Atoms { get; private set; }
        public ResidueTable Residues { get; private set; }
        public FragmentTable Fragments { get; private set; }
        public ChainTable Chains { get; private set; }
        public EntityTable Entities { get; private set; }

        ///<summary> Not to be called directly! Use ModelCollection.GetModel() or similar.</summary>
        internal Model(CifCategory cifCategory, int[] rows, int modelNumber){
            this.ModelNumber = modelNumber;
            throw new NotImplementedException();
        }
    }
}