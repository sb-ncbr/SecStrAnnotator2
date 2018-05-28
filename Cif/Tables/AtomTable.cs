using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    public class AtomTable
    {
        public const string KEY_COLUMN = ID_COLUMN;

        public readonly int Count;
        
        // down
        public readonly int[] RowIndex;

        // up
        public readonly int[] ResidueIndex;
        public readonly int[] FragmentIndex;
        public readonly int[] ChainIndex;
        public readonly int[] EntityIndex;
        private readonly Model model;

        // own properties
        public const string ID_COLUMN = "id";
        public readonly string[] Id;

        public const string NAME_COLUMN = "label_atom_id";
        public readonly string[] Name;

        public const string ELEMENT_COLUMN = "type_symbol";
        public readonly string[] Element;

        public const string ALT_LOC_COLUMN = "label_alt_id";
        public const string DEFAULT_ALT_LOC = ".";
        public readonly string[] AltLoc;

        public string String(int iAtom) => 
            model.Chains.Id[ChainIndex[iAtom]] 
            + " " + model.Residues.Compound[ResidueIndex[iAtom]] 
            + " " + model.Residues.SeqNumber[ResidueIndex[iAtom]]
            + " " + Name[iAtom]
            + (AltLoc[iAtom] != DEFAULT_ALT_LOC ? " (alt " + AltLoc[iAtom] + ")" : "" ); 

        ///<summary> Not to be called directly! Use Model.Atoms.</summary>
        internal AtomTable(Model model, CifCategory category, int[] rows, int[] atomStartsOfResidues, int[] atomStartsOfFragments, int[] atomStartsOfChains, int[] atomStartsOfEntities) {
            this.Count = rows.Length;
            // down
            this.RowIndex = rows;
            // up
            this.ResidueIndex = Model.GetUpRefs(atomStartsOfResidues);
            this.FragmentIndex = Model.GetUpRefs(atomStartsOfFragments);
            this.ChainIndex = Model.GetUpRefs(atomStartsOfChains);
            this.EntityIndex = Model.GetUpRefs(atomStartsOfEntities);
            this.model = model;
            // own properties
            RowIndex = rows;
            Id = category[ID_COLUMN].GetStrings();
            Name = category[NAME_COLUMN].GetStrings();
            Element = category[ELEMENT_COLUMN].GetStrings();
            AltLoc = category[ALT_LOC_COLUMN].GetStrings();
            // array segments
            // --
        }
    }
}