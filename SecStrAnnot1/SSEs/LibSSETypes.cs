using System;
using System.Collections.Generic;
using System.Linq;

using Cif.Components;
using protein.Libraries;

namespace protein.SSEs
{
    static class LibSSETypes
    {
        private static readonly HashSet<SSEType> ALL_HELIX_TYPES = new HashSet<SSEType>{
            SSEType.HELIX_H_TYPE,
            SSEType.HELIX_G_TYPE,
            SSEType.HELIX_I_TYPE,
            SSEType.MIXED_HELIX_TYPE
        };

        private static readonly HashSet<SSEType> ALL_SHEET_TYPES = new HashSet<SSEType>{ 
            SSEType.SHEET_TYPE, 
            SSEType.ISOLATED_BETA_BRIDGE_TYPE 
        };

        private static readonly Dictionary<SSEType, char> sseTypeToChar = new Dictionary<SSEType, char>{
            { SSEType.MIXED_HELIX_TYPE, 'h' },
            { SSEType.HELIX_G_TYPE, 'G' },
            { SSEType.HELIX_H_TYPE, 'H' },
            { SSEType.HELIX_I_TYPE, 'I' },
            { SSEType.SHEET_TYPE, 'E' },
            { SSEType.ISOLATED_BETA_BRIDGE_TYPE, 'B' },
            { SSEType.TURN_C7_TYPE, 'U' },
            { SSEType.WIGGLE_C7_TYPE, 'W' },
            { SSEType.BULGE_CLASSIC_SHORT_SIDE_TYPE, 'n' },
            { SSEType.BULGE_CLASSIC_LONG_SIDE_TYPE, 'N' },
            { SSEType.BULGE_WIDE_SHORT_SIDE_TYPE, 'm' },
            { SSEType.BULGE_WIDE_LONG_SIDE_TYPE, 'M' },
            { SSEType.BULGE_ANTIPARALLEL33_SHORT_SIDE_TYPE, 'u' }, // in 2qad chain B ~ resi 15 // "short" side is the one donating protons
            { SSEType.BULGE_ANTIPARALLEL33_LONG_SIDE_TYPE, 'U' }, // "long" side is the one accepting protons
            { SSEType.BULGE_ANTIPARALLEL22_SHORT_SIDE_TYPE, 't' }, // in 1gei ~ resi 13 // "short" side is the one donating protons
            { SSEType.BULGE_ANTIPARALLEL22_LONG_SIDE_TYPE, 'T' }, // "long" side is the one accepting protons
            { SSEType.BULGE_ANTIPARALLEL15_SHORT_SIDE_TYPE, 's' }, // in 1gjm ~ resi 94
            { SSEType.BULGE_ANTIPARALLEL15_LONG_SIDE_TYPE, 'S' },
            { SSEType.BULGE_ANTIPARALLEL23_SHORT_SIDE_TYPE, 'o' }, // in 3dbg ~ resi 452
            { SSEType.BULGE_ANTIPARALLEL23_LONG_SIDE_TYPE, 'O' },
            { SSEType.BULGE_PARALLEL14_SHORT_SIDE_TYPE, 'p' },
            { SSEType.BULGE_PARALLEL14_LONG_SIDE_TYPE, 'P' },
            { SSEType.BULGE_PARALLEL32_SHORT_SIDE_TYPE, 'q' }, // in 3ruk ~ resi 38
            { SSEType.BULGE_PARALLEL32_LONG_SIDE_TYPE, 'Q' },
            { SSEType.BULGE_PARALLEL13_SHORT_SIDE_TYPE, 'r' }, // in 3dax ~ resi 35
            { SSEType.BULGE_PARALLEL13_LONG_SIDE_TYPE, 'R' },
            { SSEType.BULGE_PARALLEL33_SHORT_SIDE_TYPE, 'l' }, // in 3v8d ~ resi 69
            { SSEType.BULGE_PARALLEL33_LONG_SIDE_TYPE, 'L' },
        };


        public static char AsChar(this SSEType type) => sseTypeToChar[type];

        public static bool IsHelix(this SSEType type) => ALL_HELIX_TYPES.Contains(type);

        public static bool IsSheet(this SSEType type) => ALL_SHEET_TYPES.Contains(type);
    }
}
