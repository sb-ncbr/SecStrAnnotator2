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
        private int[] atomStartIndex;

        private int[] residueStartIndex;

        private int[] fragmentStartIndex;

        private int[] chainStartIndex;
        public int ChainStartIndex(int iEntity) => chainStartIndex[iEntity];
        public int ChainEndIndex(int iEntity) => chainStartIndex[iEntity+1];

        // up
        // nothing

        // own properties
        private const string ID_COLUMN = "label_entity_id";
        private string[] id;
        public string Id(int iEntity) => id[iEntity];

        ///<summary> Not to be called directly! Use Model.Entities.</summary>
        public EntityTable(CifCategory category, int[] rows, out int[] entityStarts) {
            throw new NotImplementedException();
        }
    }
}