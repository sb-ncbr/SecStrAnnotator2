using System;
using System.Collections.Generic;
using System.Linq;

namespace SecStrAnnot2.Cif.Tables
{
    // Based on http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/atom_site.html
    public class AtomTable
    {
        public bool UsingAuth { get; set; } // switches between auth_* items and label_* items

        public bool[] IsHetatm { get; private set; } // [group_PDB]
        public string[] Id { get; private set; } // [id] unique atom identifier (mandatory, unique)
        public string[] ElementSymbol { get; private set; } // [type_symbol] element symbol (mandatory) => _atom_type.symbol
        private string[] labelName; // [label_atom_id] atom name, e.g. CA, ND1... (mandatory) => _chem_comp_atom.atom_id
        private string[] AltLoc; // [label_alt_id] alternative location (mandatory) => _atom_sites_alt.id
        private string[] labelResName; // [label_comp_id] residue name (mandatory) => _chem_comp.id
        private string[] labelChain; // [label_asym_id] chain ID (mandatory) => _struct_asym.id
        public string[] Entity { get; private set; } // [label_entity_id] entity ID (mandatory) => _entity.id
        private int[] labelResSeq; // [label_seq_id] residue sequence number (mandatory) => _entity_poly_seq.num
        public string[] InsertionCode; // [pdbx_PDB_ins_code] residue insertion code (mandatory) => _atom_site.label_asym_id
        public double[] X; // [Cartn_x] x coordinate
        public double[] Y; // [Cartn_y] y coordinate
        public double[] Z; // [Cartn_z] z coordinate
        public double[] Occupancy; // [occupancy] occupancy (default 1.0)
        public double[] BIso; // [B_iso_or_equiv] isotropic atomic displacement parameter, or equivalent
        public int[] Charge; // [pdbx_formal_charge] atom charge (default 0)
        private string[] authName; // [auth_atom_id] atom name, e.g. CA, ND1...
        private string[] authResName; // [auth_comp_id] residu name
        private string[] authResSeq; // [auth_seq_id] residue sequence number (string!)
        private string[] authChain; // [auth_asym_id] chain ID
        public int[] ModelNum; // [pdbx_PDB_model_num] model number
        
        public AtomTable (CifBlock cifBlock) : this(cifBlock["_atom_site"]) { }
        public AtomTable (CifCategory cifCategory){
            IsHetatm = cifCategory["group_PDB"].GetStrings().Select(s => s == "HETATM").ToArray(); //TODO catch exceptions on this and other non-mandatory items
            Id = cifCategory["id"].GetStrings();
            ElementSymbol = cifCategory["type_symbol"].GetStrings();
            labelName = cifCategory["label_atom_id"].GetStrings();
            AltLoc = cifCategory["label_alt_id"].GetStrings();
            labelResName = cifCategory["label_comp_id"].GetStrings();
            labelChain = cifCategory["label_asym_id"].GetStrings();
            Entity = cifCategory["label_entity_id"].GetStrings();
            //labelResSeq = cifCategory["label_seq_id"].GetIntegers(); //TODO catch . on hetatms (make method GetIntegersOrDefault()?)
            InsertionCode = cifCategory["pdbx_PDB_ins_code"].GetStrings();
            X = cifCategory["Cartn_x"].GetDoubles();
            Y = cifCategory["Cartn_y"].GetDoubles();
            Z = cifCategory["Cartn_z"].GetDoubles();
            Occupancy = cifCategory["occupancy"].GetDoubles();
            BIso = cifCategory["B_iso_or_equiv"].GetDoubles();
            //Charge = cifCategory["pdbx_formal_charge"].GetIntegers(); //TODO catch ? on hetatms (make method GetIntegersOrDefault()?)
            authName = cifCategory["auth_atom_id"].GetStrings();
            authResName = cifCategory["auth_comp_id"].GetStrings();
            authResSeq = cifCategory["auth_seq_id"].GetStrings();
            authChain = cifCategory["auth_asym_id"].GetStrings();
            ModelNum = cifCategory["pdbx_PDB_model_num"].GetIntegers();
        }
    }
}