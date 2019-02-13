using System;
using System.Collections.Generic;
using System.Linq;

namespace /*SecStrAnnot2.*/Cif.Tables
{
    public class FragmentTable
    {
        //TODO make arrays private and access individual values through methods
        public readonly int Count;

        // down
        private readonly int[] atomStartIndex;
        public readonly ArraySegment<int> AtomStartIndex;
        public readonly ArraySegment<int> AtomEndIndex;

        private int[] residueStartIndex;
        public readonly ArraySegment<int> ResidueStartIndex;
        public readonly ArraySegment<int> ResidueEndIndex;
        
        // up
        public readonly int[] ChainIndex;        
        public readonly int[] EntityIndex;
        private readonly Model model;

        public string String(int iFragment) => 
            model.Chains.Id[ChainIndex[iFragment]] 
            + " " + model.Residues.SeqNumber[ResidueStartIndex[iFragment]]
            + "-" + model.Residues.SeqNumber[ResidueEndIndex[iFragment] - 1]; 
        
        // own properties
        // --  


        ///<summary> Not to be called directly! Use Model.Fragments.</summary>
        internal FragmentTable(Model model, CifCategory category, int[] rows, int[] atomStartsOfFragments, int[] residueStartsOfFragments, int[] fragmentStartsOfChains, int[] fragmentStartsOfEntities) {
            this.Count = atomStartsOfFragments.Length - 1;
            // down
            this.atomStartIndex = atomStartsOfFragments;
            this.residueStartIndex = residueStartsOfFragments;
            // up
            this.ChainIndex = Model.GetUpRefs(fragmentStartsOfChains);
            this.EntityIndex = Model.GetUpRefs(fragmentStartsOfEntities);
            this.model = model;
            // own properties
            // --
            // array segments
            AtomStartIndex = new ArraySegment<int>(atomStartIndex, 0, Count);
            AtomEndIndex = new ArraySegment<int>(atomStartIndex, 1, Count);
            ResidueStartIndex = new ArraySegment<int>(residueStartIndex, 0, Count);
            ResidueEndIndex = new ArraySegment<int>(residueStartIndex, 1, Count);
        }

        ///<summary> Not to be called directly! Use Model.Fragments.</summary>
        internal FragmentTable(Model model, int[] atomStartsOfFragments, int[] residueStartsOfFragments, int[] fragmentStartsOfChains, int[] fragmentStartsOfEntities) {
            this.Count = atomStartsOfFragments.Length - 1;
            // down
            this.atomStartIndex = atomStartsOfFragments;
            this.residueStartIndex = residueStartsOfFragments;
            // up
            this.ChainIndex = Model.GetUpRefs(fragmentStartsOfChains);
            this.EntityIndex = Model.GetUpRefs(fragmentStartsOfEntities);
            this.model = model;
            // own properties
            // --
            // array segments
            AtomStartIndex = new ArraySegment<int>(atomStartIndex, 0, Count);
            AtomEndIndex = new ArraySegment<int>(atomStartIndex, 1, Count);
            ResidueStartIndex = new ArraySegment<int>(residueStartIndex, 0, Count);
            ResidueEndIndex = new ArraySegment<int>(residueStartIndex, 1, Count);
        }
    }
}