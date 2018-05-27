using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    public class ChainTable
    {
        public const string KEY_COLUMN = ID_COLUMN;
        
        //TODO make arrays private and access individual values through methods
        public int Count { get; private set; }

        // down
        private int[] atomStartIndex;

        private int[] residueStartIndex;

        private int[] fragmentStartIndex;
        public int FragmentStartIndex(int iChain) => fragmentStartIndex[iChain];
        public int FragmentEndIndex(int iChain) => fragmentStartIndex[iChain+1];
        
        // up
        private int[] entityIndex;
        public int EntityIndex(int iChain) => entityIndex[iChain];

        // own properties
        private const string ID_COLUMN = "label_asym_id";
        private string[] id;
        public string Id(int iChain) => id[iChain];

        private const string AUTH_ID_COLUMN = "auth_asym_id";
        private string[]  authId;
        public string AuthId(int iChain) => authId[iChain];
    }
}