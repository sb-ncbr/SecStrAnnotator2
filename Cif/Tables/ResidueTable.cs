using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    public class ResidueTable
    {
        public const string KEY_COLUMN = SEQ_NUMBER_COLUMN;
        public const int DEFAULT_RESIDUE_NUMBER = 0; // used for residues with label_seq_id = .

        //TODO make arrays private and access individual values through methods
        public readonly int Count;

        // down
        private readonly int[] atomStartIndex;
        public readonly ArraySegment<int> AtomStartIndex;
        public readonly ArraySegment<int> AtomEndIndex;

        // up
        public readonly int[] FragmentIndex;
        public readonly int[] ChainIndex;
        public readonly int[] EntityIndex;
        private readonly Model model;

        // own properties
        public const string SEQ_NUMBER_COLUMN = "label_seq_id";
        public readonly int[] SeqNumber;

        public const string COMPOUND_COLUMN = "label_comp_id";
        public readonly string[] Compound;

        public string String(int iResidue) => 
            model.Chains.Id[ChainIndex[iResidue]] 
            + " " + Compound[iResidue] 
            + " " + SeqNumber[iResidue]; 


        ///<summary> Not to be called directly! Use Model.Residues.</summary>
        internal ResidueTable(Model model, CifCategory category, int[] rows, int[] atomStartsOfResidues, int[] residueStartsOfFragments, int[] residueStartsOfChains, int[] residueStartsOfEntities) {
            this.Count = atomStartsOfResidues.Length - 1;
            // down
            this.atomStartIndex = atomStartsOfResidues;
            // up
            this.FragmentIndex = Model.GetUpRefs(residueStartsOfFragments);
            this.ChainIndex = Model.GetUpRefs(residueStartsOfChains);
            this.EntityIndex = Model.GetUpRefs(residueStartsOfEntities);
            this.model = model;
            // own properties
            this.SeqNumber = category[SEQ_NUMBER_COLUMN].GetIntegers(Model.GetSelectedElements(rows, atomStartsOfResidues, true), DEFAULT_RESIDUE_NUMBER);
            this.Compound = category[COMPOUND_COLUMN].GetStrings(Model.GetSelectedElements(rows, atomStartsOfResidues, true));
            // array segments
            AtomStartIndex = new ArraySegment<int>(atomStartIndex, 0, Count);
            AtomEndIndex = new ArraySegment<int>(atomStartIndex, 1, Count);
        }

        // usefull info about residues in CIF is in _pdbx_poly_seq_scheme        
    }
}