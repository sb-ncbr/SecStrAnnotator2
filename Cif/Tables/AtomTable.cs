using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    public class AtomTable
    {
        private const string KEY_COLUMN = ID_COLUMN;

        //TODO make arrays private and access individual values through methods
        public int Count { get; private set; }
        
        // down
        private int[] rowIndex;
        public int RowIndex(int iAtom) => rowIndex[iAtom];

        // up
        private int[] residueIndex;
        public int ResidueIndex(int iAtom) => residueIndex[iAtom];

        private int[] fragmentIndex;
        public int FragmentIndex(int iAtom) => fragmentIndex[iAtom];

        private int[] chainIndex;
        public int ChainIndex(int iAtom) => chainIndex[iAtom];

        private int[] entityIndex;
        public int EntityIndex(int iAtom) => entityIndex[iAtom];

        // own properties
        private const string ID_COLUMN = "id";
        private string[] id;
        public string Id(int iAtom) => id[iAtom];

        private const string NAME_COLUMN = "label_atom_id";
        private string[] name;
        public string Name(int iAtom) => name[iAtom];

        private const string ELEMENT_COLUMN = "type_symbol";
        private string[] element;
        public string Element(int iAtom) => element[iAtom];

        private const string ALT_LOC_COLUMN = "label_alt_id";
        private string[] altLoc;
        public string AltLoc(int iAtom) => altLoc[iAtom];

        internal AtomTable(CifCategory cifCategory, int[] rows, int[] residueStarts){
            rowIndex = rows;
            id = cifCategory[ID_COLUMN].GetStrings();
            name = cifCategory[NAME_COLUMN].GetStrings();
            element = cifCategory[ELEMENT_COLUMN].GetStrings();
            altLoc = cifCategory[ALT_LOC_COLUMN].GetStrings();
            //TODO connections to residues etc.
        }
    }
}