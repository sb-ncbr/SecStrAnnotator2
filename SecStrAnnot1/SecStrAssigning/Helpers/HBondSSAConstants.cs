
namespace protein.SecStrAssigning.Helpers
{
    static class HBondSSAConstants {
        public const int MIN_HBONDS_PER_LADDER = 2;
        public const int MIN_HBONDS_PER_HELIX = 2; //1;
        public const int MIN_OVERLAP_FOR_JOINING = 0; // condition for joining 2 beta-strands: end0 >= start1 + MIN_OVERLAP et vice versa
        public const bool STRANDS_BY_ALPHA = true; // if true then residues are assigned to a strand if they have C-alpha included in a cycle(DSSP-style), if false then residues are assigned to a strand if they have any atom included in a cycle(HERA-style)
        public const bool BULGES_BY_ALPHA = false;
        public const bool HELICES_BY_ALPHA = true;
        // To use DSSP-style joining, set MIN_Z_OVERLAP_FOR_JOINING = 3*MIN_OVERLAP_FOR_JOINING+2
        public const int MIN_Z_OVERLAP_FOR_JOINING = 0; // 3*MIN_OVERLAP_FOR_JOINING+2; // condition for joining 2 beta-strands: Z(end0) >= Z(start1) + MIN_OVERLAP et vice versa
        public const bool ALLOW_BULGE_A33 = true;  // antiparallel beta-bulge defined as only 1 missing H-bond from regular beta-ladder(in 2qad chain B ~ resi 15)
        public const int MAX_Z_SHIFT_ON_SHORT_BULGE_SIDE = 4;
        public const int MAX_Z_SHIFT_ON_LONG_BULGE_SIDE = 13;
    }
}