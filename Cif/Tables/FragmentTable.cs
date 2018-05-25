using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    public class FragmentTable
    {
        //TODO make arrays private and access individual values through methods
        public int Count { get; private set; }

        // down
        private int[] atomStartIndex;

        private int[] residueStartIndex;
        public int ResidueStartIndex(int iFragment) => residueStartIndex[iFragment];
        public int ResiduEndIndex(int iFragment) => residueStartIndex[iFragment+1];
        
        // up
        private int[] chainIndex;
        public int ChainIndex(int iFragment) => chainIndex[iFragment];

        private int[] entityIndex;
        public int EntityIndex(int iFragment) => entityIndex[iFragment];

        // own properties
        // --  
    }
}