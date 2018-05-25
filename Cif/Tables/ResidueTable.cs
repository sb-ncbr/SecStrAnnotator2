using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    public class ResidueTable
    {
        private const string KEY_COLUMN = SEQ_NUMBER_COLUMN;

        //TODO make arrays private and access individual values through methods
        public int Count { get; private set; }

        // down
        private int[] atomStartIndex;
        public int AtomStartIndex(int iResidue) => atomStartIndex[iResidue];
        public int AtomEndIndex(int iResidue) => atomStartIndex[iResidue+1];

        // up
        private int[] fragmentIndex;
        public int FragmentIndex(int iResidue) => fragmentIndex[iResidue];

        private int[] chainIndex;
        public int ChainIndex(int iResidue) => chainIndex[iResidue];

        private int[] entityIndex;
        public int EntityIndex(int iResidue) => entityIndex[iResidue];

        // own properties
        private const string SEQ_NUMBER_COLUMN = "label_seq_id";
        private int[] seqNumber;
        public int SeqNumber(int iResidue) => seqNumber[iResidue];

        private const string COMPOUND_COLUMN = "label_comp_id";
        private string[] compound;
        public string Compound(int iResidue) => compound[iResidue];

        // usefull info about residues in CIF is in _pdbx_poly_seq_scheme        
    }
}