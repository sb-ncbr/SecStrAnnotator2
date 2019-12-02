
namespace protein.SecStrAssigning.Helpers
{
    enum MicrostrandType{		
        REGULAR_PARALEL,  // 'E';
        REGULAR_ANTIPARALEL,  // 'E';

		BULGE_UNSPECIFIED_PARALLEL_SHORT_SIDE,
		BULGE_UNSPECIFIED_PARALLEL_LONG_SIDE,
		BULGE_UNSPECIFIED_PARALLEL_EQUAL_SIDE,
		BULGE_UNSPECIFIED_ANTIPARALLEL_SHORT_SIDE,
		BULGE_UNSPECIFIED_ANTIPARALLEL_LONG_SIDE,
		BULGE_UNSPECIFIED_ANTIPARALLEL_EQUAL_SIDE,

		BULGE_CLASSIC_SHORT_SIDE,  // 'n';
		BULGE_CLASSIC_LONG_SIDE,  // 'N';
		BULGE_WIDE_SHORT_SIDE,  // 'm';
		BULGE_WIDE_LONG_SIDE,  // 'M';
		BULGE_ANTIPARALLEL33_SHORT_SIDE,  // 'u'; // in 2qad chain B ~ resi 15 // "short" side is the one donating protons
		BULGE_ANTIPARALLEL33_LONG_SIDE,  // 'U'; // "long" side is the one accepting protons
		BULGE_ANTIPARALLEL22_SHORT_SIDE,  // 't'; // in 1gei ~ resi 13 // "short" side is the one donating protons
		BULGE_ANTIPARALLEL22_LONG_SIDE,  // 'T'; // "long" side is the one accepting protons
		BULGE_ANTIPARALLEL15_SHORT_SIDE,  // 's'; // in 1gjm ~ resi 94
		BULGE_ANTIPARALLEL15_LONG_SIDE,  // 'S'; 
		BULGE_ANTIPARALLEL23_SHORT_SIDE,  // 'o'; // in 3dbg ~ resi 452
		BULGE_ANTIPARALLEL23_LONG_SIDE,  // 'O'; 
		BULGE_PARALLEL14_SHORT_SIDE,  // 'p';
		BULGE_PARALLEL14_LONG_SIDE,  // 'P';
		BULGE_PARALLEL32_SHORT_SIDE,  // 'q'; // in 3ruk ~ resi 38
		BULGE_PARALLEL32_LONG_SIDE,  // 'Q';
		BULGE_PARALLEL13_SHORT_SIDE,  // 'r'; // in 3dax ~ resi 35
		BULGE_PARALLEL13_LONG_SIDE,  // 'R';
		BULGE_PARALLEL33_SHORT_SIDE,  // 'l'; // in 3v8d ~ resi 69
		BULGE_PARALLEL33_LONG_SIDE,  // 'L';
    }
}