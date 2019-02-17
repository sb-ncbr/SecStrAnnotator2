using System;
using System.Collections.Generic;
using System.Linq;

namespace /*SecStrAnnot2.*/Cif.Tables
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

        // own optional properties
        public const string AUTH_SEQ_NUMBER_COLUMN = "auth_seq_id";
        public readonly MaybeValueArray<int> AuthSeqNumber;

        public const string AUTH_INSERTION_CODE_COLUMN = "pdbx_PDB_ins_code";
        public readonly MaybeClassArray<string> AuthInsertionCode;

        public const string AUTH_COMPOUND_COLUMN = "auth_comp_id";
        public readonly MaybeClassArray<string> AuthCompound;

        public string String(int iResidue) => 
            model.Chains.Id[ChainIndex[iResidue]] 
            + " " + Compound[iResidue] 
            + " " + SeqNumber[iResidue]; 


        ///<summary> Not to be called directly! Use Model.Residues.</summary>
        internal ResidueTable(Model model, CifCategory category, int[] rows, 
                              int[] atomStartsOfResidues, int[] residueStartsOfFragments, int[] residueStartsOfChains, int[] residueStartsOfEntities) {
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
            // own optional properties
            this.AuthSeqNumber = new MaybeValueArray<int>(category.ContainsItem(AUTH_SEQ_NUMBER_COLUMN) ? category[AUTH_SEQ_NUMBER_COLUMN].GetIntegers(Model.GetSelectedElements(rows, atomStartsOfResidues, true), DEFAULT_RESIDUE_NUMBER) : null);
            this.AuthInsertionCode = new MaybeClassArray<string>(category.ContainsItem(AUTH_INSERTION_CODE_COLUMN) ? category[AUTH_INSERTION_CODE_COLUMN].GetStrings(Model.GetSelectedElements(rows, atomStartsOfResidues, true)) : null);
            this.AuthCompound = new MaybeClassArray<string>(category.ContainsItem(AUTH_COMPOUND_COLUMN) ? category[AUTH_COMPOUND_COLUMN].GetStrings(Model.GetSelectedElements(rows, atomStartsOfResidues, true)) : null);
            // array segments
            AtomStartIndex = new ArraySegment<int>(atomStartIndex, 0, Count);
            AtomEndIndex = new ArraySegment<int>(atomStartIndex, 1, Count);
        }

        ///<summary> Not to be called directly! Use Model.Residues.</summary>
        internal ResidueTable(Model model, int[] atomStartsOfResidues, int[] residueStartsOfFragments, int[] residueStartsOfChains, int[] residueStartsOfEntities,
                              ResidueInfo[] residueInfo){
            this.Count = atomStartsOfResidues.Length - 1;
            // down
            this.atomStartIndex = atomStartsOfResidues;
            // up
            this.FragmentIndex = Model.GetUpRefs(residueStartsOfFragments);
            this.ChainIndex = Model.GetUpRefs(residueStartsOfChains);
            this.EntityIndex = Model.GetUpRefs(residueStartsOfEntities);
            this.model = model;
            // own properties
            this.SeqNumber = residueInfo.Select(i => i.SeqNumber).ToArray();
            this.Compound = residueInfo.Select(i => i.Compound).ToArray();
            this.AuthSeqNumber = new MaybeValueArray<int>(residueInfo.Select(i => i.AuthSeqNumber).ToArray());
            this.AuthInsertionCode = new MaybeClassArray<string>(residueInfo.Select(i => i.AuthInsertionCode).ToArray());
            this.AuthCompound = new MaybeClassArray<string>(residueInfo.Select(i => i.AuthCompound).ToArray());
            // own optional properties
            this.AuthSeqNumber = new MaybeValueArray<int>(residueInfo.Select(i => i.AuthSeqNumber).ToArray());
            this.AuthInsertionCode = new MaybeClassArray<string>(residueInfo.Select(i => i.AuthInsertionCode).ToArray());
            this.AuthCompound = new MaybeClassArray<string>(residueInfo.Select(i => i.AuthCompound).ToArray());
            // array segments
            AtomStartIndex = new ArraySegment<int>(atomStartIndex, 0, Count);
            AtomEndIndex = new ArraySegment<int>(atomStartIndex, 1, Count);
        }

        // usefull info about residues in CIF is in _pdbx_poly_seq_scheme        
    }
}