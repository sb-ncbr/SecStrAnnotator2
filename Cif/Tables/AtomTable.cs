using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    // Based on http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/atom_site.html
    public class AtomTable
    {
        private const string CATEGORY_NAME = "_atom_site";
        public const int DEFAULT_RES_SEQ = 0; // assign label_seq_id 0 to residues without label_seq_id (mmCIF allows only positive values of label_seq_id)
        public const int DEFAULT_CHARGE = 0;

        public int Count { get; private set; } // number of atoms

        // Mandatory data items
        public string[] Id { get; private set; } // [id] unique atom identifier (mandatory, unique)
        public string[] authChain { get; private set; } // [auth_asym_id] chain ID
        public string[] AltLoc { get; private set; } // [label_alt_id] alternative location (mandatory) => _atom_sites_alt.id
        public string[] labelChain { get; private set; } // [label_asym_id] chain ID (mandatory) => _struct_asym.id
        public string[] labelName { get; private set; } // [label_atom_id] atom name, e.g. CA, ND1... (mandatory) => _chem_comp_atom.atom_id
        public string[] labelResName { get; private set; } // [label_comp_id] residue name (mandatory) => _chem_comp.id
        public string[] Entity { get; private set; } // [label_entity_id] entity ID (mandatory) => _entity.id
        public int[] labelResSeq { get; private set; } // [label_seq_id] residue sequence number (mandatory) => _entity_poly_seq.num
        public string[] ElementSymbol { get; private set; } // [type_symbol] element symbol (mandatory) => _atom_type.symbol
        // Optional data items
        public string[] authName { get; private set; } // [auth_atom_id] atom name, e.g. CA, ND1...
        public string[] authResName { get; private set; } // [auth_comp_id] residu name
        public string[] authResSeq_String { get; private set; } // [auth_seq_id] residue sequence number (string!)
        public double[] BIso { get; private set; } // [B_iso_or_equiv] isotropic atomic displacement parameter, or equivalent
        public double[] X { get; private set; } // [Cartn_x] x coordinate
        public double[] Y { get; private set; } // [Cartn_y] y coordinate
        public double[] Z { get; private set; } // [Cartn_z] z coordinate
        public bool[] IsHetatm { get; private set; } // [group_PDB]
        public double[] Occupancy { get; private set; } // [occupancy] occupancy (default 1.0)
        public int[] Charge { get; private set; } // [pdbx_formal_charge] atom charge (default 0)
        public string[] InsertionCode { get; private set; } // [pdbx_PDB_ins_code] residue insertion code (mandatory) => _atom_site.label_asym_id
        public int[] ModelNum { get; private set; } // [pdbx_PDB_model_num] model number
        // Processed items
        private int[] _authResSeq_Int; // [auth_seq_id] parsed to int
        public int[] authResSeq_Int { 
            get { 
                if (_authResSeq_Int == null) {
                    TryParseAuthResSeq(); 
                } 
                return _authResSeq_Int; 
            } 
        } // [auth_seq_id] parsed to int


        /// <summary> Switches between auth_* items and label_* items (*atom_id, *seq_id, *comp_id, *asym_id).</summary>
        public bool UsingAuth { get; set; }

        public string[] Name => UsingAuth ? authName : labelName;
        public string[] ResName => UsingAuth ? authResName : labelResName;
        public int[] ResSeq => UsingAuth ? authResSeq_Int : labelResSeq;
        public string[] Chain => UsingAuth ? authChain : labelChain;
        

        public AtomTable (CifBlock cifBlock) : this(cifBlock[CATEGORY_NAME]) { }

        public AtomTable (CifBlock cifBlock, int[] rows) : this(cifBlock[CATEGORY_NAME], rows) { }

        public AtomTable (CifCategory cifCategory){
            // Mandatory data items
            try {
                Id = cifCategory["id"].GetStrings();
                authChain = cifCategory["auth_asym_id"].GetStrings();
                AltLoc = cifCategory["label_alt_id"].GetStrings();
                labelChain = cifCategory["label_asym_id"].GetStrings();
                labelName = cifCategory["label_atom_id"].GetStrings();
                labelResName = cifCategory["label_comp_id"].GetStrings();
                Entity = cifCategory["label_entity_id"].GetStrings();
                labelResSeq = cifCategory["label_seq_id"].GetIntegers(DEFAULT_RES_SEQ);
                ElementSymbol = cifCategory["type_symbol"].GetStrings();
            } catch (KeyNotFoundExceptionWithKey<string> e){
                Lib.WriteError("Mandatory data item " + CATEGORY_NAME + "." + e.Key + " is missing in the CIF file");
                throw e;
            }
            // Optional data items
            try {
                authName = cifCategory["auth_atom_id"].GetStrings();
                authResName = cifCategory["auth_comp_id"].GetStrings();
                authResSeq_String = cifCategory["auth_seq_id"].GetStrings();
                BIso = cifCategory["B_iso_or_equiv"].GetDoubles();
                X = cifCategory["Cartn_x"].GetDoubles();
                Y = cifCategory["Cartn_y"].GetDoubles();
                Z = cifCategory["Cartn_z"].GetDoubles();
                IsHetatm = cifCategory["group_PDB"].GetStrings().Select(s => s == "HETATM").ToArray();
                Occupancy = cifCategory["occupancy"].GetDoubles();
                Charge = cifCategory["pdbx_formal_charge"].GetIntegers(DEFAULT_CHARGE);
                InsertionCode = cifCategory["pdbx_PDB_ins_code"].GetStrings();
                ModelNum = cifCategory["pdbx_PDB_model_num"].GetIntegers();
            } catch (KeyNotFoundExceptionWithKey<string> e){
                Lib.WriteWarning(CATEGORY_NAME + "." + e.Key + " is missing in the CIF file");
            }
            // Atom count
            Count = Id.Length;
        }

        public AtomTable (CifCategory cifCategory, int[] rows){
            // Mandatory data items
            try {
                Id = cifCategory["id"].GetStrings(rows);
                authChain = cifCategory["auth_asym_id"].GetStrings(rows);
                AltLoc = cifCategory["label_alt_id"].GetStrings(rows);
                labelChain = cifCategory["label_asym_id"].GetStrings(rows);
                labelName = cifCategory["label_atom_id"].GetStrings(rows);
                labelResName = cifCategory["label_comp_id"].GetStrings(rows);
                Entity = cifCategory["label_entity_id"].GetStrings(rows);
                labelResSeq = cifCategory["label_seq_id"].GetIntegers(rows, DEFAULT_RES_SEQ);
                ElementSymbol = cifCategory["type_symbol"].GetStrings(rows);
            } catch (KeyNotFoundExceptionWithKey<string> e){
                Lib.WriteError("Mandatory data item " + CATEGORY_NAME + "." + e.Key + " is missing in the CIF file");
                throw e;
            }
            // Optional data items
            try {
                authName = cifCategory["auth_atom_id"].GetStrings(rows);
                authResName = cifCategory["auth_comp_id"].GetStrings(rows);
                authResSeq_String = cifCategory["auth_seq_id"].GetStrings(rows);
                BIso = cifCategory["B_iso_or_equiv"].GetDoubles(rows);
                X = cifCategory["Cartn_x"].GetDoubles(rows);
                Y = cifCategory["Cartn_y"].GetDoubles(rows);
                Z = cifCategory["Cartn_z"].GetDoubles(rows);
                IsHetatm = cifCategory["group_PDB"].GetStrings(rows).Select(s => s == "HETATM").ToArray();
                Occupancy = cifCategory["occupancy"].GetDoubles(rows);
                Charge = cifCategory["pdbx_formal_charge"].GetIntegers(rows, DEFAULT_CHARGE);
                InsertionCode = cifCategory["pdbx_PDB_ins_code"].GetStrings(rows);
                ModelNum = cifCategory["pdbx_PDB_model_num"].GetIntegers(rows);
            } catch (KeyNotFoundExceptionWithKey<string> e){
                Lib.WriteWarning(CATEGORY_NAME + "." + e.Key + " is missing in the CIF file");
            }
            // Atom count
            Count = Id.Length;
        }

        private void TryParseAuthResSeq(){
            int n = authResSeq_String.Length;
            int[] ints = new int[n];
            for (int i = 0; i < n; i++)
            {
                bool ok = int.TryParse(authResSeq_String[i], out ints[i]);
                if (!ok) {
                    if (authResSeq_String[i] == "." || authResSeq_String[i] == "?"){
                        ints[i] = DEFAULT_RES_SEQ;
                    }
                }
            }
        }


    }
}