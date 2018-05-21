namespace SecStrAnnot2.Cif.Tables
{
    public class ResidueTable
    {
        // usefull info about residues in CIF is in _pdbx_poly_seq_scheme

        public int Count { get; private set; } // number of atoms

        // Mandatory data items
        public int[] Id { get; private set; } // [_atom_site.label_seq_id]
        public string[] authReqSeq { get; private set; } // [_atom_site.label_seq_id]
        
    }
}