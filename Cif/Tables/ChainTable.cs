using System;
using System.Collections.Generic;
using System.Linq;

namespace /*SecStrAnnot2.*/Cif.Tables
{
    public class ChainTable
    {
        public const string KEY_COLUMN = ID_COLUMN;
        
        //TODO make arrays private and access individual values through methods
        public readonly int Count;

        // down
        private readonly int[] atomStartIndex;
        public readonly ArraySegment<int> AtomStartIndex;
        public readonly ArraySegment<int> AtomEndIndex;

        private int[] residueStartIndex;
        public readonly ArraySegment<int> ResidueStartIndex;
        public readonly ArraySegment<int> ResidueEndIndex;

        private int[] fragmentStartIndex;
        public readonly ArraySegment<int> FragmentStartIndex;
        public readonly ArraySegment<int> FragmentEndIndex;

        // up
        public readonly int[] EntityIndex;
        private readonly Model model;

        // own properties
        public const string ID_COLUMN = "label_asym_id";
        public readonly string[] Id;

        public const string AUTH_ID_COLUMN = "auth_asym_id";
        public readonly string[] AuthId;

        public string String(int iChain) => Id[iChain];


        ///<summary> Not to be called directly! Use Model.Chains.</summary>
        internal ChainTable(Model model, CifCategory category, int[] rows, int[] atomStartsOfChains, int[] residueStartsOfChains, int[] fragmentStartsOfChains, int[] chainStartsOfEntities) {
            this.Count = atomStartsOfChains.Length - 1;
            // down
            this.atomStartIndex = atomStartsOfChains;
            this.residueStartIndex = residueStartsOfChains;
            this.fragmentStartIndex = fragmentStartsOfChains;
            // up
            this.EntityIndex = Model.GetUpRefs(chainStartsOfEntities);
            this.model = model;
            // own properties
            this.Id = category[ID_COLUMN].GetStrings(Model.GetSelectedElements(rows, atomStartsOfChains, true));
            this.AuthId = category[AUTH_ID_COLUMN].GetStrings(Model.GetSelectedElements(rows, atomStartsOfChains, true));
            // array segments
            AtomStartIndex = new ArraySegment<int>(atomStartIndex, 0, Count);
            AtomEndIndex = new ArraySegment<int>(atomStartIndex, 1, Count);
            ResidueStartIndex = new ArraySegment<int>(residueStartIndex, 0, Count);
            ResidueEndIndex = new ArraySegment<int>(residueStartIndex, 1, Count);
            FragmentStartIndex = new ArraySegment<int>(fragmentStartIndex, 0, Count);
            FragmentEndIndex = new ArraySegment<int>(fragmentStartIndex, 1, Count);
        }

        ///<summary> Not to be called directly! Use Model.Chains.</summary>
        internal ChainTable(Model model, int[] atomStartsOfChains, int[] residueStartsOfChains, int[] fragmentStartsOfChains, int[] chainStartsOfEntities,
                            string[] id, string[] authId) {
            this.Count = atomStartsOfChains.Length - 1;
            // down
            this.atomStartIndex = atomStartsOfChains;
            this.residueStartIndex = residueStartsOfChains;
            this.fragmentStartIndex = fragmentStartsOfChains;
            // up
            this.EntityIndex = Model.GetUpRefs(chainStartsOfEntities);
            this.model = model;
            // own properties
            this.Id = id;
            this.AuthId = authId;
            // array segments
            AtomStartIndex = new ArraySegment<int>(atomStartIndex, 0, Count);
            AtomEndIndex = new ArraySegment<int>(atomStartIndex, 1, Count);
            ResidueStartIndex = new ArraySegment<int>(residueStartIndex, 0, Count);
            ResidueEndIndex = new ArraySegment<int>(residueStartIndex, 1, Count);
            FragmentStartIndex = new ArraySegment<int>(fragmentStartIndex, 0, Count);
            FragmentEndIndex = new ArraySegment<int>(fragmentStartIndex, 1, Count);
        }
    }
}