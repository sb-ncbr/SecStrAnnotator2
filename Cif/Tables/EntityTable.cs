using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    public class EntityTable
    {
        public const string KEY_COLUMN = ID_COLUMN;

        //TODO make arrays private and access individual values through methods
        public int Count { get; private set; }

        // down
        internal int[] atomStartIndex;

        internal int[] residueStartIndex;

        internal int[] fragmentStartIndex;

        internal int[] chainStartIndex;
        public int ChainStartIndex(int iEntity) => chainStartIndex[iEntity];
        public int ChainEndIndex(int iEntity) => chainStartIndex[iEntity+1];

        // up
        // --

        // own properties
        private const string ID_COLUMN = "label_entity_id";
        private string[] id;
        public string Id(int iEntity) => id[iEntity];

        ///<summary> Not to be called directly! Use Model.Entities.</summary>
        public EntityTable(CifCategory category, int[] rows, int[] atomStartsOfEntities, int[] residueStartsOfEntities, int[] fragmentStartsOfEntities, int[] chainStartsOfEntities) {
            this.Count = atomStartsOfEntities.Length - 1;
            // down
            this.atomStartIndex = atomStartsOfEntities;
            this.residueStartIndex = residueStartsOfEntities;
            this.fragmentStartIndex = fragmentStartsOfEntities;
            this.chainStartIndex = chainStartsOfEntities;
            // up
            // --
            // own properties
            this.id = category[KEY_COLUMN].GetStrings(Model.GetSelectedElements(rows, atomStartsOfEntities, true));
        }
    }
}