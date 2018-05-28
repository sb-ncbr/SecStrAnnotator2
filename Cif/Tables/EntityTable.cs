using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    public class EntityTable
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

        private int[] chainStartIndex;
        public readonly ArraySegment<int> ChainStartIndex;
        public readonly ArraySegment<int> ChainEndIndex;
        // public int ChainStartIndex(int iEntity) => chainStartIndex[iEntity];
        // public int ChainEndIndex(int iEntity) => chainStartIndex[iEntity+1];

        // up
        private readonly Model model;

        // own properties
        public const string ID_COLUMN = "label_entity_id";
        public readonly string[] Id;

        public string String(int iEntity) => Id[iEntity];


        ///<summary> Not to be called directly! Use Model.Entities.</summary>
        internal EntityTable(Model model, CifCategory category, int[] rows, int[] atomStartsOfEntities, int[] residueStartsOfEntities, int[] fragmentStartsOfEntities, int[] chainStartsOfEntities) {
            Count = atomStartsOfEntities.Length - 1;
            // down
            atomStartIndex = atomStartsOfEntities;
            residueStartIndex = residueStartsOfEntities;
            fragmentStartIndex = fragmentStartsOfEntities;
            chainStartIndex = chainStartsOfEntities;
            // up
            this.model = model;
            // own properties
            Id = category[ID_COLUMN].GetStrings(Model.GetSelectedElements(rows, atomStartsOfEntities, true));
            // array segments
            AtomStartIndex = new ArraySegment<int>(atomStartIndex, 0, Count);
            AtomEndIndex = new ArraySegment<int>(atomStartIndex, 1, Count);
            ResidueStartIndex = new ArraySegment<int>(residueStartIndex, 0, Count);
            ResidueEndIndex = new ArraySegment<int>(residueStartIndex, 1, Count);
            FragmentStartIndex = new ArraySegment<int>(fragmentStartIndex, 0, Count);
            FragmentEndIndex = new ArraySegment<int>(fragmentStartIndex, 1, Count);
            ChainStartIndex = new ArraySegment<int>(chainStartIndex, 0, Count);
            ChainEndIndex = new ArraySegment<int>(chainStartIndex, 1, Count);
        }
    }
}