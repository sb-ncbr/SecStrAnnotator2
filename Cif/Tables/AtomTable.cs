using System;
using System.Collections.Generic;
using System.Linq;
using Cif.Libraries;

namespace Cif.Tables
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

        public const string IS_HETATM_COLUMN = "group_PDB";
        public readonly bool[] IsHetatm;

        public const string X_COLUMN = "Cartn_x";
        public readonly double[] X;
        public const string Y_COLUMN = "Cartn_y";
        public readonly double[] Y;
        public const string Z_COLUMN = "Cartn_z";
        public readonly double[] Z;

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
            Id = category[ID_COLUMN].GetStrings(rows);
            Name = category[NAME_COLUMN].GetStrings(rows);
            Element = category[ELEMENT_COLUMN].GetStrings(rows);
            AltLoc = category[ALT_LOC_COLUMN].GetStrings(rows);
            try { X = category[X_COLUMN].GetDoubles(rows); } catch (KeyNotFoundExceptionWithKey<string> e) { Lib.WriteWarning("Missing data item " + e.Key); }
            try { Y = category[Y_COLUMN].GetDoubles(rows); } catch (KeyNotFoundExceptionWithKey<string> e) { Lib.WriteWarning("Missing data item " + e.Key); }
            try { Z = category[Z_COLUMN].GetDoubles(rows); } catch (KeyNotFoundExceptionWithKey<string> e) { Lib.WriteWarning("Missing data item " + e.Key); }
            try { IsHetatm = category[IS_HETATM_COLUMN].GetTrueWhereFirstCharacterMatches(rows, 'H'); } catch (KeyNotFoundExceptionWithKey<string> e) { Lib.WriteWarning("Missing data item " + e.Key); }
            // array segments
            // --
        }

        ///<summary> Not to be called directly! Use Model.Atoms.</summary>
        internal AtomTable(Model model, int[] atomStartsOfResidues, int[] atomStartsOfFragments, int[] atomStartsOfChains, int[] atomStartsOfEntities, 
                           string[] id, AtomInfo[] info) {
            this.Count = info.Length;
            // down
            this.RowIndex = Enumerable.Range(0, Count).ToArray();
            // up
            this.ResidueIndex = Model.GetUpRefs(atomStartsOfResidues);
            this.FragmentIndex = Model.GetUpRefs(atomStartsOfFragments);
            this.ChainIndex = Model.GetUpRefs(atomStartsOfChains);
            this.EntityIndex = Model.GetUpRefs(atomStartsOfEntities);
            this.model = model;
            // own properties
            Id = id;
            Name = info.Select(i => i.Name).ToArray();
            Element = info.Select(i => i.Element).ToArray();
            AltLoc = info.Select(i => i.AltLoc).ToArray();
            IsHetatm = info.Select(i => i.IsHetatm).ToArray();
            X = info.Select(i => i.X).ToArray();
            Y = info.Select(i => i.Y).ToArray();
            Z = info.Select(i => i.Z).ToArray();
            // array segments
            // --
        }
    }
}